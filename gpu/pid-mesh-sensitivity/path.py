#!/usr/bin/env python3
"""Lateral path of the bubble in the lab frame (path instability).

Usage: path.py [--csv summary.csv] <run dir>...

The PID controller keeps the bubble at the centre of the mesh by moving the
frame, so a sideways drift of the bubble shows up in the lab-frame velocity
rise_u, rise_w of data.csv (each row the mean over its interval) rather than
in its position. The lab-frame lateral displacement is their integral from
t = 0, in bubble diameters; like drift.py in gpu/mesh-sensitivity, this
reports when it first exceeds 0.01, 0.05 and 0.25, the exponential growth
rate of the lateral speed fitted over 1e-4 < speed < 1e-2 (where it grows
from noise), the speed at t = 0 that fit extrapolates to, and the drift
direction at the end of the fit.
"""
import argparse
import csv
import math
import os

import numpy as np

THRESHOLDS = (0.01, 0.05, 0.25)
FIT_RANGE = (1e-4, 1e-2)


def load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def lateral_path(run_dir):
    """Time, lab-frame lateral displacement (x, z) and lateral speed."""
    d = load(os.path.join(run_dir, "data.csv"))
    t = d["Time"]
    dt = np.diff(t, prepend=0.0)
    x = np.cumsum(d["rise_u"]*dt)
    z = np.cumsum(d["rise_w"]*dt)
    return t, x, z, np.hypot(d["rise_u"], d["rise_w"])


def summarize(run_dir):
    t, x, z, speed = lateral_path(run_dir)
    r = np.hypot(x, z)
    out = {"run": os.path.basename(os.path.normpath(run_dir)), "t_end": t[-1], "r_end": r[-1],
           "speed_end": speed[-1]}
    for thr in THRESHOLDS:
        above = np.nonzero(r > thr)[0]
        out["t_r>{:g}".format(thr)] = t[above[0]] if len(above) else float("nan")

    # First contiguous stretch of samples inside the fit range, after the
    # start-up transient (t > 1).
    inside = (speed > FIT_RANGE[0]) & (speed < FIT_RANGE[1]) & (t > 1.0)
    idx = np.nonzero(inside)[0]
    if len(idx) > 2:
        stop = idx[0] + np.argmax(~inside[idx[0]:]) if (~inside[idx[0]:]).any() else len(t)
        sel = np.arange(idx[0], stop)
    else:
        sel = np.array([], dtype=int)
    if len(sel) > 2:
        rate, log_s0 = np.polyfit(t[sel], np.log(speed[sel]), 1)
        out.update(growth_rate=rate, fit_t0=t[sel[0]], fit_t1=t[sel[-1]], speed0=math.exp(log_s0),
                   direction_deg=math.degrees(math.atan2(z[sel[-1]], x[sel[-1]])))
    else:
        out.update(growth_rate=float("nan"), fit_t0=float("nan"), fit_t1=float("nan"),
                   speed0=float("nan"), direction_deg=float("nan"))
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--csv", help="also write the summary to this CSV file")
    parser.add_argument("dirs", nargs="+")
    args = parser.parse_args()

    results = [summarize(d) for d in args.dirs]
    keys = ["t_r>{:g}".format(thr) for thr in THRESHOLDS]
    print("{:>8} {:>6} {:>7} {:>8} {:>8} {:>8} {:>8} {:>7} {:>11} {:>8} {:>5}".format(
        "run", "t_end", "r_end", "speed", *keys, "rate", "fit t", "speed(0)", "dir"))
    for r in results:
        print("{:>8} {:6.2f} {:7.3f} {:8.1e} {:8.2f} {:8.2f} {:8.2f} {:7.3f} {:5.1f}-{:<5.1f} {:8.1e} {:5.0f}".format(
            r["run"], r["t_end"], r["r_end"], r["speed_end"], *[r[k] for k in keys], r["growth_rate"],
            r["fit_t0"], r["fit_t1"], r["speed0"], r["direction_deg"]))
    if args.csv:
        with open(args.csv, "w") as f:
            w = csv.DictWriter(f, fieldnames=list(results[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(results)


if __name__ == "__main__":
    main()
