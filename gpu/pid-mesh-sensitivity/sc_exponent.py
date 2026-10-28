#!/usr/bin/env python3
"""Apparent Schmidt number exponent n in Sh ~ Sc^n from two species of a run.

Usage: sc_exponent.py [--tmin 10] [--tmax 30] [--species c4] [--trace 0.5] <run dir>...

c (data.csv, at the Sc of case.hpp) and a further species c<Sc>
(data-c<Sc>.csv) see the same flow, so the ratio of their Sh cancels most of
what the flow does to both, e.g. a bubble's shape oscillation. n is
ln(Sh_c<Sc>/Sh_c)/ln(Sc/Sc_c) of the window means of Sh (rows weighted by
their interval, as in sherwood.py), and +/- the standard error of its 1 time
unit batch values. n_sink is the same from Sh_sink. --trace also prints n
over consecutive intervals of that length from the start of the species, to
see when one started afresh (speciesStart) has settled.
"""
import argparse
import csv
import math
import os
import re

import numpy as np


def load(path):
    with open(path) as f:
        rows = list(csv.DictReader(f))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def case_sc(run_dir):
    with open(os.path.join(run_dir, "case.hpp")) as f:
        m = re.search(r"^static double Sc = ([0-9.eE+-]+);", f.read(), re.M)
    if not m:
        raise SystemExit("{}: no Sc in case.hpp".format(run_dir))
    return float(m.group(1))


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--tmin", type=float, default=10.0)
    parser.add_argument("--tmax", type=float, default=30.0)
    parser.add_argument("--species", default="c4", help="the further species c<Sc> (default c4)")
    parser.add_argument("--trace", type=float, default=0.0,
                        help="also print n over intervals of this length (default: no trace)")
    parser.add_argument("dirs", nargs="+")
    args = parser.parse_args()
    sc = float(args.species[1:])

    lines = []
    print("{:>8} {:>6} {:>9} {:>7} {:>7} {:>6} {:>6} {:>7}".format(
        "run", "t_end", "samples", "Sh_c", "Sh_" + args.species, "n", "+/-", "n_sink"))
    for d in args.dirs:
        log_sc = math.log(sc/case_sc(d))
        c = load(os.path.join(d, "data.csv"))
        s = load(os.path.join(d, "data-{}.csv".format(args.species)))
        # data.csv may start earlier (with the row a run started by setup-runs.sh
        # FROM copies), so take its rows at the further species' times (the
        # same step writes the same time into both files).
        t = s["Time"]
        rows = np.minimum(np.searchsorted(c["Time"], t - 1e-9), len(c["Time"]) - 1)
        if not np.allclose(c["Time"][rows], t, rtol=0, atol=1e-9):
            raise SystemExit("{}: data.csv lacks rows of data-{}.csv".format(d, args.species))
        weight = np.diff(t, prepend=np.nan)
        interval = np.nanmedian(weight)
        # The species' first row ends the step after its start, data.csv's
        # previous row.
        weight[0] = t[0] - c["Time"][rows[0] - 1] if rows[0] > 0 else interval

        def n_of(a, b, sh="Sh"):
            """n, the number of rows and the two mean Sh over (a, b]."""
            inside = (t > a + 0.5*interval) & (t <= b + 0.5*interval)
            if not inside.any():
                return float("nan"), 0, float("nan"), float("nan")
            w = weight[inside]
            mean_c = (w*c[sh][rows][inside]).sum()/w.sum()
            mean_s = (w*s[sh][inside]).sum()/w.sum()
            return math.log(mean_s/mean_c)/log_sc, int(inside.sum()), mean_c, mean_s

        n, samples, sh_c, sh_s = n_of(args.tmin, args.tmax)
        # Complete batches only.
        batches = [n_of(b, b + 1) for b in np.arange(args.tmin, args.tmax - 1e-9)]
        batches = np.array([x[0] for x in batches if x[1] == int(round(1/interval))])
        se = batches.std(ddof=1)/math.sqrt(len(batches)) if len(batches) > 1 else float("nan")
        print("{:>8} {:6.2f} {:4d}/{:<4d} {:7.3f} {:7.3f} {:6.3f} {:6.3f} {:7.3f}".format(
            os.path.basename(os.path.normpath(d)), t[-1], samples,
            int(round((args.tmax - args.tmin)/interval)), sh_c, sh_s, n, se,
            n_of(args.tmin, args.tmax, "Sh_sink")[0]))
        if args.trace > 0:
            start = t[0] - weight[0]
            k = 0
            while start + (k + 1)*args.trace <= t[-1] + 0.5*interval:
                a = start + k*args.trace
                lines.append("{:>8} {:6.2f}-{:<6.2f} {:6.3f}".format(
                    os.path.basename(os.path.normpath(d)), a, a + args.trace,
                    n_of(a, a + args.trace)[0]))
                k += 1
    if lines:
        print()
        print("{:>8} {:>13} {:>6}".format("run", "t", "n"))
        print("\n".join(lines))


if __name__ == "__main__":
    main()
