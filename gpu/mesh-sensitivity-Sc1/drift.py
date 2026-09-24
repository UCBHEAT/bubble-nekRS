#!/usr/bin/env python3
"""Lateral drift of the bubble away from its initial rise axis.

Usage: drift.py [--csv summary.csv] <run dir>...

The bubble centroid comes from <run dir>/data.csv (x_c and z_c, written by
bubble3d.udf every sampleInterval in double precision) or, for runs without
those columns, <run dir>/post.csv (animate.py, from single precision
checkpoints, so offsets below about 1e-6 are noise). The lateral offset r
from the initial axis, in bubble diameters, grows exponentially once the
bubble starts to leave its periodic wake; this reports when r first exceeds
0.01, 0.05 and 0.25, the growth rate fitted over 1e-3 < r < 0.1, the offset
at t = 0 that fit extrapolates to, and the drift direction.
"""
import argparse
import csv
import math
import os

import numpy as np

L = 4.0  # periodic domain width in x and z (bubble3d.box)
THRESHOLDS = (0.01, 0.05, 0.25)
FIT_RANGE = (1e-3, 0.1)


def load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def unwrap(x):
    return L/(2*math.pi)*np.unwrap(2*math.pi*x/L)


def centroid(run_dir):
    path = os.path.join(run_dir, "data.csv")
    if os.path.exists(path):
        d = load(path)
        if "x_c" in d:
            return "data.csv", d["Time"], unwrap(d["x_c"]), unwrap(d["z_c"])
    d = load(os.path.join(run_dir, "post.csv"))
    return "post.csv", d["Time"], d["x_c_unwrapped"], d["z_c_unwrapped"]


def summarize(run_dir):
    source, t, x, z = centroid(run_dir)
    dx, dz = x - x[0], z - z[0]
    r = np.hypot(dx, dz)
    out = {"run": os.path.basename(os.path.normpath(run_dir)), "source": source,
           "t_end": t[-1], "r_end": r[-1]}
    for thr in THRESHOLDS:
        above = np.nonzero(r > thr)[0]
        out["t_r>{:g}".format(thr)] = t[above[0]] if len(above) else float("nan")

    # First contiguous stretch of samples inside the fit range.
    inside = (r > FIT_RANGE[0]) & (r < FIT_RANGE[1])
    idx = np.nonzero(inside)[0]
    if len(idx) > 2:
        stop = idx[0] + np.argmax(~inside[idx[0]:]) if (~inside[idx[0]:]).any() else len(r)
        sel = np.arange(idx[0], stop)
        rate, log_r0 = np.polyfit(t[sel], np.log(r[sel]), 1)
        out.update(growth_rate=rate, fit_t0=t[sel[0]], fit_t1=t[sel[-1]], r0=math.exp(log_r0),
                   direction_deg=math.degrees(math.atan2(dz[sel[-1]], dx[sel[-1]])))
    else:
        out.update(growth_rate=float("nan"), fit_t0=float("nan"), fit_t1=float("nan"),
                   r0=float("nan"), direction_deg=float("nan"))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--csv", help="also write the summary to this CSV file")
    parser.add_argument("dirs", nargs="+")
    args = parser.parse_args()

    results = [summarize(d) for d in args.dirs]
    keys = ["t_r>{:g}".format(thr) for thr in THRESHOLDS]
    print("{:>8} {:>8} {:>6} {:>6} {:>8} {:>8} {:>8} {:>7} {:>11} {:>8} {:>5}".format(
        "run", "source", "t_end", "r_end", *keys, "rate", "fit t", "r(0)", "dir"))
    for r in results:
        print("{:>8} {:>8} {:6.2f} {:6.3f} {:8.2f} {:8.2f} {:8.2f} {:7.3f} {:5.1f}-{:<5.1f} {:8.1e} {:5.0f}".format(
            r["run"], r["source"], r["t_end"], r["r_end"], *[r[k] for k in keys], r["growth_rate"],
            r["fit_t0"], r["fit_t1"], r["r0"], r["direction_deg"]))
    if args.csv:
        with open(args.csv, "w") as f:
            w = csv.DictWriter(f, fieldnames=list(results[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(results)


if __name__ == "__main__":
    main()
