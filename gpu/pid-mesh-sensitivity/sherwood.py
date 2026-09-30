#!/usr/bin/env python3
"""Average the Sherwood number and the bubble's rise over a time window for
each run directory.

Usage: sherwood.py [--tmin 10] [--tmax 30] [--csv summary.csv] <run dir>...

Each data.csv row holds mdot and the rise velocity averaged over the sample
interval ending at Time, with total_area sampled at Time, so the window
average uses the rows whose intervals lie inside (tmin, tmax]. Rows are
weighted by the length of their interval, so the short rows at job
boundaries (bubble.udf also samples each job's last step) count
proportionally less. The uncertainty is the standard error of 1 time unit
batch means, which are much less correlated than individual samples.

Sh is the species budget (mdot, see bubble.udf); Sh_sink is what the hard
sink counted, the way gpu/mesh-sensitivity measured it. Re_rise is Re times
the lab-frame rise speed and the volume-equivalent diameter, the Reynolds
number the correlations need; u_lateral is the mean lab-frame lateral speed
of the bubble (0 while it rises straight).
"""
import argparse
import csv
import os
import re

import numpy as np


def case_numbers(run_dir):
    """Re and Sc from the run's case.hpp."""
    with open(os.path.join(run_dir, "case.hpp")) as f:
        text = f.read()
    value = {}
    for name in ("Re", "Sc"):
        m = re.search(r"^static double {} = ([0-9.eE+-]+);".format(name), text, re.M)
        if not m:
            raise SystemExit("{}: no {} in case.hpp".format(run_dir, name))
        value[name] = float(m.group(1))
    return value["Re"], value["Sc"]


def load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def summarize(run_dir, tmin, tmax):
    Re, Sc = case_numbers(run_dir)
    data = load(os.path.join(run_dir, "data.csv"))
    t = data["Time"]
    weight = np.diff(t, prepend=0.0)
    interval = np.median(weight)
    sel = (t > tmin + 0.5*interval) & (t <= tmax + 0.5*interval)
    n_expected = int(round((tmax - tmin)/interval))
    w = weight[sel]

    def mean(x):
        return (w*x[sel]).sum()/w.sum()

    def batch_se(x):
        means = []
        for b in np.arange(tmin, tmax):
            inside = (t[sel] > b + 0.5*interval) & (t[sel] <= b + 1 + 0.5*interval)
            if inside.any():
                means.append((w[inside]*x[sel][inside]).sum()/w[inside].sum())
        means = np.array(means)
        return means.std(ddof=1)/np.sqrt(len(means)) if len(means) > 1 else float("nan")

    sh = data["Sh"]
    sh_mean = mean(sh)
    u_rise = mean(np.hypot(np.hypot(data["rise_u"], data["rise_v"]), data["rise_w"]))
    d_eq = (6.0*mean(data["gas_volume"])/np.pi)**(1.0/3.0)
    return {
        "run": os.path.basename(os.path.normpath(run_dir)),
        "t_end": t[-1],
        "samples": int(sel.sum()),
        "expected_samples": n_expected,
        "Sh_mean": sh_mean,
        "Sh_std": np.sqrt((w*(sh[sel] - sh_mean)**2).sum()/w.sum()) if sel.sum() > 1 else float("nan"),
        "Sh_batch_se": batch_se(sh),
        "Sh_ratio_of_means": Re*Sc*mean(data["mdot"])/mean(data["total_area"]),
        "Sh_sink_mean": mean(data["Sh_sink"]),
        "area_mean": mean(data["total_area"]),
        "c_bulk_mean": mean(data["c_bulk"]),
        "rise_velocity": mean(data["rise_v"]),
        "u_lateral": mean(np.hypot(data["rise_u"], data["rise_w"])),
        "d_eq": d_eq,
        "Re_rise": Re*u_rise*d_eq,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tmin", type=float, default=10.0)
    parser.add_argument("--tmax", type=float, default=30.0)
    parser.add_argument("--csv", help="also write the summary to this CSV file")
    parser.add_argument("dirs", nargs="+")
    args = parser.parse_args()

    results = [summarize(d, args.tmin, args.tmax) for d in args.dirs]
    print("{:>8} {:>6} {:>9} {:>7} {:>6} {:>6} {:>9} {:>7} {:>6} {:>6} {:>7} {:>7} {:>5}".format(
        "run", "t_end", "samples", "Sh", "+/-", "std", "Sh(means)", "Sh_sink", "area", "c_bulk",
        "rise", "lateral", "Re"))
    for r in results:
        print("{run:>8} {t_end:6.2f} {samples:4d}/{expected_samples:<4d} {Sh_mean:7.3f} "
              "{Sh_batch_se:6.3f} {Sh_std:6.3f} {Sh_ratio_of_means:9.3f} {Sh_sink_mean:7.3f} "
              "{area_mean:6.3f} {c_bulk_mean:6.4f} {rise_velocity:7.4f} {u_lateral:7.4f} "
              "{Re_rise:5.0f}".format(**r))
    if args.csv:
        with open(args.csv, "w") as f:
            w = csv.DictWriter(f, fieldnames=list(results[0]), lineterminator="\n")
            w.writeheader()
            w.writerows(results)


if __name__ == "__main__":
    main()
