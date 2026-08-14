#!/usr/bin/env python3
"""Build an HDF5 PS time series + regression velocity GeoTIFF from GMTSAR
displacement outputs, for any region (generalized from the Kaohsiung-only
make_pt_ts_h5.py so the same script works for other counties/AOIs).

Steps:
  1. Read per-epoch displacement files disp_NNN_ll.xy (lon lat disp_mm)
     and epoch list data_date.txt (YYYYMMDD, one per disp file) from
     --ts-dir.
  2. Keep only PS points inside the --boundary county/AOI polygon.
  2b. If --vel-geojson is given, re-trend each PS series to the calibrated
     velocity field (Point features, field_3 = mm/yr): the raw disp series
     may carry an uncorrected orbital/reference ramp, so per point
     D += (v_cal - v_ols) * t.  After this the OLS regression of every PS
     series equals its calibrated velocity.  Points without a coordinate
     match in the velocity file use the plane fit of (v_ols - v_cal) as
     the correction.  --deramp optionally removes a regional plane from
     the calibrated velocity field itself before it is used (see below).
     If --vel-geojson is not given, this step is skipped entirely and a
     warning is printed.
  3. Moving-window ordinary kriging of each epoch onto a regular
     lon/lat grid: each grid cell is solved from its K nearest PS
     points (batched LAPACK solves).  pykrige is used only to fit the
     per-epoch variogram on a random subsample.  PS coordinates are
     identical for every epoch, so neighbor indices and distance
     matrices are computed once.
  4. Clip the grid with the boundary polygon (outside -> NaN); only
     in-polygon cells are solved.  Cells farther than --max-ps-dist-km
     from the nearest PS point are also NaN, so data gaps (mountains)
     are left blank instead of being extrapolated (set to 0 to disable).
  5. Write an HDF5 file with the same schema as
     levelingmap/output_data/leveling_ts_2021_2026.h5:
       /date        (n,)  bytes string
       /timeseries  (n, rows, cols) float32, NaN outside mask
       root attrs: UNIT, DISPLAY_UNIT, EPSG, GEOTRANSFORM
  6. Write the per-pixel OLS regression velocity of the final grid
     series as a GeoTIFF (same grid) for gmtsar_timeseries_html_viewer.py --vel,
     so the mapped velocity is mathematically identical to the
     regression the viewer computes when a cell is clicked.

Example (Kaohsiung, matches the original make_pt_ts_h5.py output):
  python3 gmtsar2h5.py \
    --ts-dir 202504-202605 --boundary KS.geojson \
    --vel-geojson ks202504-202604.geojson \
    --deramp -18.1124,48.6986,-74.9069,120.3629025741,22.7490741316 \
    --out-h5 output_data/KS_202504-202605TS.h5 \
    --out-vel output_data/KS_202504-202605_vel.tif

Run with the isce2 conda env (needs pykrige, shapely, scipy, h5py, gdal).
"""

import os

# Batched LAPACK solves over many small (K+1)x(K+1) matrices spin-lock and
# contend for cores under multi-threaded BLAS; force single-threaded BLAS
# before numpy (which loads the BLAS backend) is imported.
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import argparse
import json
import sys
import time
from datetime import datetime

import h5py
import numpy as np
import shapely
from osgeo import gdal, osr
from pykrige.ok import OrdinaryKriging
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist
from shapely.geometry import shape
from shapely.ops import unary_union

KM_PER_DEG = 111.19    # local equirectangular projection scale
VARIO_SAMPLE = 5000    # subsample size for per-epoch variogram fitting
CELL_BATCH = 100000    # grid cells per batched LAPACK solve


def parse_args():
    p = argparse.ArgumentParser(
        description="Build an HDF5 PS time series + regression velocity "
                    "GeoTIFF from GMTSAR disp_NNN_ll.xy outputs, clipped "
                    "to an AOI polygon, with optional re-trend to a "
                    "calibrated velocity field.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    p.add_argument("--ts-dir", required=True,
                    help="folder with disp_NNN_ll.xy and data_date.txt")
    p.add_argument("--boundary", required=False, default=None,
                    help="AOI polygon geojson (county/city boundary); "
                         "未給時不裁切 (all PS points treated as inside)")
    p.add_argument("--out-h5", required=True, help="output HDF5 path")
    p.add_argument("--out-vel", required=True,
                    help="output regression-velocity GeoTIFF path")
    p.add_argument("--vel-geojson", default=None,
                    help="calibrated velocity field, Point features with "
                         "field_3 = mm/yr; if omitted the retrend step is "
                         "skipped and a warning is printed")
    p.add_argument("--deramp", default=None,
                    help="C0,C1,C2,LON0,LAT0: subtract plane "
                         "C0+C1*(lon-LON0)+C2*(lat-LAT0) (mm/yr, "
                         "mm/yr/deg) from --vel-geojson before use; only "
                         "effective together with --vel-geojson")
    p.add_argument("--grid-step", type=float, default=0.001,
                    help="output grid spacing, deg")
    p.add_argument("--margin", type=float, default=0.005,
                    help="padding around PS extent, deg")
    p.add_argument("--max-ps-dist-km", type=float, default=0.5,
                    help="grid cells farther than this from the nearest "
                         "PS point are left NaN (data-gap mask); 0 "
                         "disables the mask")
    p.add_argument("--k-neigh", type=int, default=64,
                    help="PS points per moving-window kriging solve")

    # argparse only recognizes "-18.1" as a negative number (not an
    # option); "--deramp -18.1,...,-74.9,..." has commas, so it fails
    # that check and argparse treats the value as an unrecognized
    # option.  Rewrite "--deramp VALUE" to "--deramp=VALUE" so the
    # plain space-separated CLI form still works.
    argv = sys.argv[1:]
    fixed = []
    i = 0
    while i < len(argv):
        if argv[i] == "--deramp" and i + 1 < len(argv):
            fixed.append(f"--deramp={argv[i + 1]}")
            i += 2
        else:
            fixed.append(argv[i])
            i += 1
    return p.parse_args(fixed)


def read_dates(ts_dir):
    with open(os.path.join(ts_dir, "data_date.txt")) as f:
        return [ln.strip() for ln in f if ln.strip()]


def read_disp(ts_dir, idx):
    """Return (lon, lat, disp_mm) arrays for epoch index (1-based)."""
    path = os.path.join(ts_dir, f"disp_{idx:03d}_ll.xy")
    arr = np.loadtxt(path)
    return arr[:, 0], arr[:, 1], arr[:, 2]


def load_boundary(geojson_path):
    with open(geojson_path) as f:
        gj = json.load(f)
    boundary = unary_union([shape(ft["geometry"]) for ft in gj["features"]])
    # shapely.contains_xy() is called twice on this geometry (PS points,
    # grid points); without a prepared geometry each call falls back to
    # an unindexed predicate that is ~1000x slower for the same result.
    shapely.prepare(boundary)
    return boundary


def build_grid(lon, lat, grid_step, margin):
    """Regular grid (pixel centers) covering the PS extent plus margin."""
    x0 = lon.min() - margin
    x1 = lon.max() + margin
    y0 = lat.min() - margin
    y1 = lat.max() + margin
    ncols = int(np.ceil((x1 - x0) / grid_step))
    nrows = int(np.ceil((y1 - y0) / grid_step))
    # pixel centers; row 0 = north (GDAL north-up convention)
    xc = x0 + (np.arange(ncols) + 0.5) * grid_step
    yc = y1 - (np.arange(nrows) + 0.5) * grid_step
    geotransform = np.array([x0, grid_step, 0.0, y1, 0.0, -grid_step])
    return xc, yc, geotransform


def to_km(lon, lat, lat0):
    """Local equirectangular projection (deg -> km)."""
    x = np.asarray(lon) * KM_PER_DEG * np.cos(np.radians(lat0))
    y = np.asarray(lat) * KM_PER_DEG
    return np.column_stack([x, y])


def decyears(dates):
    """Years elapsed since the first epoch (day count / 365.25)."""
    d0 = datetime.strptime(dates[0], "%Y%m%d")
    return np.array([(datetime.strptime(d, "%Y%m%d") - d0).days / 365.25
                     for d in dates])


def ols_slope(t, D):
    """OLS slope along axis 0 of D (n_epochs x n_points) vs time t (yr)."""
    tm = t - t.mean()
    return (tm[:, None] * (D - D.mean(0))).sum(0) / (tm ** 2).sum()


def retrend_to_calibrated(lon, lat, D, t, vel_gj, deramp_mm,
                          deramp_lon0, deramp_lat0):
    """Force the OLS trend of each PS series to its calibrated velocity.

    Reads vel_gj (Point features, field_3 = calibrated LOS velocity in
    mm/yr), matches points to (lon, lat) by rounded coordinates, and adds
    (v_cal - v_ols) * t to each series.  Unmatched points fall back to the
    plane fit of (v_cal - v_ols) evaluated at their coordinates.
    """
    with open(vel_gj) as f:
        gj = json.load(f)
    glon = np.array([ft["geometry"]["coordinates"][0]
                     for ft in gj["features"]])
    glat = np.array([ft["geometry"]["coordinates"][1]
                     for ft in gj["features"]])
    gv = np.array([ft["properties"]["field_3"] for ft in gj["features"]])
    gv = gv - (deramp_mm[0]
               + deramp_mm[1] * (glon - deramp_lon0)
               + deramp_mm[2] * (glat - deramp_lat0))
    print(f"calibrated velocity after deramp: p1/p50/p99 = "
          f"{np.round(np.percentile(gv, [1, 50, 99]), 1)} mm/yr", flush=True)

    v_ols = ols_slope(t, D)
    lut = {(round(x, 7), round(y, 7)): i for i, (x, y) in
           enumerate(zip(glon, glat))}
    v_cal = np.full(len(lon), np.nan)
    for i, (x, y) in enumerate(zip(lon, lat)):
        j = lut.get((round(x, 7), round(y, 7)))
        if j is not None:
            v_cal[i] = gv[j]
    matched = np.isfinite(v_cal)

    # plane fit of the correction on matched points, for the unmatched rest
    corr = v_cal - v_ols
    G = np.column_stack([np.ones(matched.sum()),
                         lon[matched] - lon.mean(),
                         lat[matched] - lat.mean()])
    coef, *_ = np.linalg.lstsq(G, corr[matched], rcond=None)
    resid_std = float(np.std(corr[matched] - G @ coef))
    if not matched.all():
        Gu = np.column_stack([np.ones((~matched).sum()),
                              lon[~matched] - lon.mean(),
                              lat[~matched] - lat.mean()])
        corr[~matched] = Gu @ coef
    print(f"re-trend to calibrated velocity: matched "
          f"{matched.sum()}/{len(lon)} PS, correction plane "
          f"(const,dlon,dlat)=({coef[0]:+.1f},{coef[1]:+.1f},{coef[2]:+.1f})"
          f" mm/yr, plane residual std {resid_std:.2f} mm/yr", flush=True)
    return D + corr[None, :] * t[:, None]


def write_vel_tif(path, vel_grid, gt):
    """Write the regression velocity grid as float32 GeoTIFF (EPSG:4326)."""
    drv = gdal.GetDriverByName("GTiff")
    ds = drv.Create(path, vel_grid.shape[1], vel_grid.shape[0], 1,
                    gdal.GDT_Float32, options=["COMPRESS=DEFLATE"])
    ds.SetGeoTransform(tuple(gt))
    srs = osr.SpatialReference()
    srs.ImportFromEPSG(4326)
    ds.SetProjection(srs.ExportToWkt())
    band = ds.GetRasterBand(1)
    band.SetNoDataValue(float("nan"))
    band.WriteArray(vel_grid.astype(np.float32))
    ds = None


class MovingWindowKriging:
    """Per-cell OK from the K nearest PS points, shared geometry.

    Neighbor indices, target-neighbor distances and neighbor pairwise
    distances depend only on coordinates, so they are computed once;
    each epoch only refits the variogram and re-solves the batched
    (K+1)x(K+1) OK systems.
    """

    def __init__(self, pts_km, targets_km, k):
        self.k = k
        tree = cKDTree(pts_km)
        self.d_tgt, self.idx = tree.query(targets_km, k=k, workers=-1)
        print(f"  neighbor search done: {len(targets_km)} cells, "
              f"tgt-dist p50/p95 {np.median(self.d_tgt[:, 0]):.3f}/"
              f"{np.percentile(self.d_tgt[:, 0], 95):.3f} km", flush=True)
        m = len(targets_km)
        self.d_pair = np.empty((m, k, k), dtype=np.float32)
        for s in range(0, m, CELL_BATCH):
            e = min(s + CELL_BATCH, m)
            nb = pts_km[self.idx[s:e]]              # (b, k, 2)
            diff = nb[:, :, None, :] - nb[:, None, :, :]
            self.d_pair[s:e] = np.sqrt((diff ** 2).sum(-1))
        print("  pairwise neighbor distances cached "
              f"({self.d_pair.nbytes / 1e9:.1f} GB)", flush=True)

    def fit_variogram(self, pts_km, vals, rng):
        sub = rng.choice(len(vals), min(VARIO_SAMPLE, len(vals)),
                         replace=False)
        ok = OrdinaryKriging(
            pts_km[sub, 0], pts_km[sub, 1], vals[sub],
            variogram_model="spherical", nlags=20,
            coordinates_type="euclidean",
        )
        return ok.variogram_model_parameters, ok.variogram_function

    def solve_epoch(self, pts_km, vals, rng):
        params, vfunc = self.fit_variogram(pts_km, vals, rng)
        k = self.k
        m = len(self.idx)
        out = np.empty(m, dtype=np.float32)
        for s in range(0, m, CELL_BATCH):
            e = min(s + CELL_BATCH, m)
            b = e - s
            a = np.empty((b, k + 1, k + 1))
            a[:, :k, :k] = vfunc(params, self.d_pair[s:e])
            a[:, np.arange(k), np.arange(k)] = 0.0
            a[:, k, :] = 1.0
            a[:, :, k] = 1.0
            a[:, k, k] = 0.0
            rhs = np.empty((b, k + 1))
            rhs[:, :k] = vfunc(params, self.d_tgt[s:e])
            rhs[:, k] = 1.0
            w = np.linalg.solve(a, rhs[..., None])[:, :k, 0]
            out[s:e] = (w * vals[self.idx[s:e]]).sum(1)
        return out, params


def main():
    args = parse_args()
    t0 = time.time()

    deramp_mm = (0.0, 0.0, 0.0)
    deramp_lon0 = 0.0
    deramp_lat0 = 0.0
    if args.deramp:
        parts = [float(x) for x in args.deramp.split(",")]
        assert len(parts) == 5, \
            "--deramp needs 5 comma-separated values: C0,C1,C2,LON0,LAT0"
        deramp_mm = tuple(parts[:3])
        deramp_lon0, deramp_lat0 = parts[3], parts[4]

    dates = read_dates(args.ts_dir)
    n_files = len([f for f in os.listdir(args.ts_dir)
                   if f.startswith("disp_") and f.endswith("_ll.xy")])
    assert len(dates) == n_files, f"{len(dates)} dates vs {n_files} disp files"

    # common PS coordinates come from epoch 2 (epoch 1 is the all-zero
    # reference epoch with a larger, unfiltered point list)
    lon_all, lat_all, _ = read_disp(args.ts_dir, 2)
    if args.boundary:
        boundary = load_boundary(args.boundary)
        inside = shapely.contains_xy(boundary, lon_all, lat_all)
    else:
        boundary = None
        inside = np.ones(len(lon_all), dtype=bool)
        print("未給 --boundary: 不做行政區裁切", flush=True)
    lon, lat = lon_all[inside], lat_all[inside]
    print(f"PS points inside boundary polygon: {inside.sum()} / {len(lon_all)}",
          flush=True)

    xc, yc, gt = build_grid(lon, lat, args.grid_step, args.margin)
    xx, yy = np.meshgrid(xc, yc)
    if boundary is not None:
        mask = shapely.contains_xy(boundary, xx.ravel(), yy.ravel()) \
            .reshape(xx.shape)
    else:
        mask = np.ones(xx.shape, dtype=bool)
    print(f"grid {len(yc)} rows x {len(xc)} cols, "
          f"{mask.sum()} of {mask.size} cells inside boundary polygon",
          flush=True)

    # drop cells with no PS point within max_ps_dist_km: kriging there
    # would confidently extrapolate into data gaps (e.g. mountain areas)
    lat0 = float(np.mean(lat))
    pts_km = to_km(lon, lat, lat0)
    if args.max_ps_dist_km > 0:
        d_near, _ = cKDTree(pts_km).query(
            to_km(xx[mask], yy[mask], lat0), k=1, workers=-1)
        n_before = int(mask.sum())
        mask[mask] = d_near <= args.max_ps_dist_km
        print(f"data-gap mask (> {args.max_ps_dist_km*1000:.0f} m from "
              f"nearest PS): {n_before} -> {mask.sum()} cells", flush=True)
    rows, cols = np.nonzero(mask)

    # load every epoch (inside-polygon points only), then optionally
    # re-trend the series to the calibrated velocity field before kriging
    t_yr = decyears(dates)
    D = np.zeros((len(dates), int(inside.sum())))
    for i in range(2, len(dates) + 1):
        lo, la, v = read_disp(args.ts_dir, i)
        assert np.array_equal(lo, lon_all) and np.array_equal(la, lat_all), \
            f"epoch {i}: coordinates differ from epoch 2"
        D[i - 1] = v[inside]
    if args.vel_geojson:
        D = retrend_to_calibrated(lon, lat, D, t_yr, args.vel_geojson,
                                  deramp_mm, deramp_lon0, deramp_lat0)
    else:
        print("警告：未提供 --vel-geojson，未做速度校正 (retrend skipped)",
              flush=True)

    solver = MovingWindowKriging(pts_km, to_km(xc[cols], yc[rows], lat0),
                                 args.k_neigh)

    rng = np.random.default_rng(42)
    ts = np.full((len(dates), len(yc), len(xc)), np.nan, dtype=np.float32)
    ts[0][mask] = 0.0  # reference epoch: zero displacement by definition

    for i in range(2, len(dates) + 1):
        te = time.time()
        z, params = solver.solve_epoch(pts_km, D[i - 1], rng)
        ts[i - 1][mask] = z
        print(f"epoch {i:2d}/{len(dates)} {dates[i-1]}: "
              f"range {np.nanmin(z):+.1f} .. {np.nanmax(z):+.1f} mm, "
              f"variogram(psill,range,nugget)="
              f"({params[0]:.1f},{params[1]:.2f},{params[2]:.1f}), "
              f"{time.time()-te:.1f}s", flush=True)

    with h5py.File(args.out_h5, "w") as f:
        f.create_dataset("date", data=np.array([d.encode() for d in dates]))
        f.create_dataset("timeseries", data=ts, compression="gzip")
        f.attrs["UNIT"] = "mm"
        f.attrs["DISPLAY_UNIT"] = "mm"
        f.attrs["EPSG"] = 4326
        f.attrs["GEOTRANSFORM"] = gt
    print(f"written: {args.out_h5}", flush=True)

    # per-pixel OLS regression velocity of the final grid series -> tif
    vel_grid = np.full(mask.shape, np.nan)
    vel_grid[mask] = ols_slope(t_yr, ts[:, mask].astype(np.float64))
    write_vel_tif(args.out_vel, vel_grid, gt)
    print(f"written: {args.out_vel}, regression velocity "
          f"p1/p50/p99 = {np.round(np.nanpercentile(vel_grid, [1, 50, 99]), 1)}"
          f" mm/yr ({time.time()-t0:.0f}s total)", flush=True)


if __name__ == "__main__":
    main()
