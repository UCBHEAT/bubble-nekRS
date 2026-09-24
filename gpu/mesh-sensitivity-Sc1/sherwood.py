#!/usr/bin/env python3
"""Average the Sherwood number over a time window for each run directory.

Usage: sherwood.py [--tmin 20] [--tmax 30] [--csv summary.csv] <run dir>...

Each data.csv row holds mdot averaged over the sample interval ending at
Time (the checkpoint interval in the 2.25x, 1.5x and 1x runs), with
total_area and c_bulk sampled at Time, so the window average uses the rows
whose intervals lie inside (tmin, tmax]. Rows are weighted by the length of
their interval, so the short rows at job boundaries (bubble3d.udf also
samples each job's last step) count proportionally less. The uncertainty is
the standard error of 1 time unit batch means, which are much less correlated
than individual samples.
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
    weight = np.diff(t, prepend=0.0)
    interval = np.median(weight)
    sel = (t > tmin + 0.5*interval) & (t <= tmax + 0.5*interval)
    n_expected = int(round((tmax - tmin)/interval))
    w = weight[sel]

    def mean(x):
        return (w*x[sel]).sum()/w.sum()

    sh = data["Sh"][sel]
    batch_means = []
    for b in np.arange(tmin, tmax):
        inside = (t[sel] > b + 0.5*interval) & (t[sel] <= b + 1 + 0.5*interval)
        if inside.any():
            batch_means.append((w[inside]*sh[inside]).sum()/w[inside].sum())
    batch_means = np.array(batch_means)
    se = batch_means.std(ddof=1)/np.sqrt(len(batch_means)) if len(batch_means) > 1 else float("nan")
    sh_mean = mean(data["Sh"])
    ratio_of_means = Pe*mean(data["mdot"])/(mean(data["total_area"])*mean(data["c_bulk"]))
    return {
        "run": os.path.basename(os.path.normpath(run_dir)),
        "t_end": t[-1],
        "samples": int(sel.sum()),
        "expected_samples": n_expected,
        "Sh_mean": sh_mean,
        "Sh_std": np.sqrt((w*(sh - sh_mean)**2).sum()/w.sum()) if len(sh) > 1 else float("nan"),
        "Sh_batch_se": se,
        "Sh_ratio_of_means": ratio_of_means,
        "area_mean": mean(data["total_area"]),
        "c_bulk_mean": mean(data["c_bulk"]),
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
            w = csv.DictWriter(f, fieldnames=list(results[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(results)


if __name__ == "__main__":
    main()
