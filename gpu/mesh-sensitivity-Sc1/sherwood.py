#!/usr/bin/env python3
"""Average the Sherwood number over a time window for each run directory.

Usage: sherwood.py [--tmin 20] [--tmax 30] [--csv summary.csv] <run dir>...

Each data.csv row holds mdot averaged over the checkpoint interval ending at
Time, with total_area and c_bulk sampled at Time, so the window average uses
the rows whose intervals lie inside (tmin, tmax]. The uncertainty is the
standard error of 1 time unit batch means, which are much less correlated
than individual checkpoints.
"""
import argparse
import csv
import os

import numpy as np

Pe = 231.4  # Re*Sc, Sc = 1 (case.hpp)


def load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def summarize(run_dir, tmin, tmax):
    data = load(os.path.join(run_dir, "data.csv"))
    t = data["Time"]
    interval = np.median(np.diff(t))
    sel = (t > tmin + 0.5*interval) & (t <= tmax + 0.5*interval)
    n_expected = int(round((tmax - tmin)/interval))
    sh = data["Sh"][sel]
    batches = [sh[(t[sel] > b + 0.5*interval) & (t[sel] <= b + 1 + 0.5*interval)]
               for b in np.arange(tmin, tmax)]
    batch_means = np.array([b.mean() for b in batches if len(b)])
    se = batch_means.std(ddof=1)/np.sqrt(len(batch_means)) if len(batch_means) > 1 else float("nan")
    ratio_of_means = Pe*data["mdot"][sel].mean()/(data["total_area"][sel].mean()*data["c_bulk"][sel].mean())
    return {
        "run": os.path.basename(os.path.normpath(run_dir)),
        "t_end": t[-1],
        "samples": int(sel.sum()),
        "expected_samples": n_expected,
        "Sh_mean": sh.mean(),
        "Sh_std": sh.std(ddof=1) if len(sh) > 1 else float("nan"),
        "Sh_batch_se": se,
        "Sh_ratio_of_means": ratio_of_means,
        "area_mean": data["total_area"][sel].mean(),
        "c_bulk_mean": data["c_bulk"][sel].mean(),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tmin", type=float, default=20.0)
    parser.add_argument("--tmax", type=float, default=30.0)
    parser.add_argument("--csv", help="also write the summary to this CSV file")
    parser.add_argument("dirs", nargs="+")
    args = parser.parse_args()

    results = [summarize(d, args.tmin, args.tmax) for d in args.dirs]
    print("{:>8} {:>7} {:>9} {:>8} {:>7} {:>8} {:>9} {:>7} {:>7}".format(
        "run", "t_end", "samples", "Sh", "+/-", "Sh_std", "Sh(means)", "area", "c_bulk"))
    for r in results:
        print("{run:>8} {t_end:7.2f} {samples:4d}/{expected_samples:<4d} {Sh_mean:8.3f} "
              "{Sh_batch_se:7.3f} {Sh_std:8.3f} {Sh_ratio_of_means:9.3f} {area_mean:7.3f} "
              "{c_bulk_mean:7.4f}".format(**r))
    if args.csv:
        with open(args.csv, "w") as f:
            w = csv.DictWriter(f, fieldnames=list(results[0]))
            w.writeheader()
            w.writerows(results)


if __name__ == "__main__":
    main()
