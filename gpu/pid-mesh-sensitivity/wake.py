#!/usr/bin/env python3
"""Bubble shape and the liquid velocity on the axis behind it, from a checkpoint.

Usage: wake.py [--radius R] <checkpoint> [<checkpoint with coordinates>]

nekRS writes the mesh coordinates only into each job's first checkpoint, so
for a later checkpoint pass that one too (same job, so the element order
matches). Both files must be Nek5000 #std with FP32 data, scalars tls, cls, c.

The PID frame moves with the bubble, so the velocity in the file is the
bubble-frame velocity: the liquid enters at the top at the frame velocity
(about minus the rise velocity) and flows down past the bubble. This prints
the bubble's front and rear on the vertical axis through the gas centroid and
its equatorial radius, where psi crosses 0.5 (interpolated between the GLL
points nearest the axis or the centroid's plane), their aspect ratio, and the
profile along the axis below the bubble: the vertical velocity over the
inflow speed, and c, in 0.05 D bins, each from the points nearest the axis
(up to 1.5 times the distance of the nearest, or within --radius). Liquid
moving up towards the bubble there is a standing eddy, the one part of the
wake whose concentration would renew by diffusion rather than by the flow,
which matters for how long Sh takes to settle at high Sc.
"""
import argparse

import numpy as np


def read(path):
    """Element ids, fields and time of a Nek5000 #std FP32 file."""
    with open(path, "rb") as f:
        words = f.read(132).decode().split()
        marker = np.frombuffer(f.read(4), dtype="<f4")[0]
    if abs(marker - 6.54321) > 1e-5:
        raise SystemExit(f"{path}: not little-endian FP32 #std")
    nxyz = int(words[2])*int(words[3])*int(words[4])
    nel, time, code = int(words[5]), float(words[7]), words[11]
    ids = np.fromfile(path, dtype="<i4", count=nel, offset=136)
    out, offset, i = {}, 136 + 4*nel, 0
    while i < len(code):
        if code[i] in "XU":
            # Vectors are stored element by element: x, y, z of each element.
            out[code[i]] = np.fromfile(path, dtype="<f4", count=3*nel*nxyz,
                                       offset=offset).reshape(nel, 3, nxyz)
            offset += 12*nel*nxyz
            i += 1
        elif code[i] == "P":
            offset += 4*nel*nxyz
            i += 1
        elif code[i] == "S":
            for k in range(int(code[i + 1:i + 3])):
                out[f"S{k + 1}"] = np.fromfile(path, dtype="<f4", count=nel*nxyz,
                                               offset=offset).reshape(nel, nxyz)
                offset += 4*nel*nxyz
            i += 3
        else:
            i += 1
    return ids, out, time


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--radius", type=float,
                        help="largest distance from the axis of the points in the profile "
                             "(default 1.5 times that of the nearest points in each bin)")
    parser.add_argument("checkpoint")
    parser.add_argument("coordinates", nargs="?", help="first checkpoint of the same job")
    args = parser.parse_args()

    ids, fields, time = read(args.checkpoint)
    if "X" not in fields:
        if not args.coordinates:
            raise SystemExit(f"{args.checkpoint} has no coordinates; pass the job's first checkpoint too")
        ids_x, fields_x, _ = read(args.coordinates)
        if not np.array_equal(ids, ids_x):
            raise SystemExit("element order differs: use the first checkpoint of the same job")
        fields["X"] = fields_x["X"]
    x, y, z = (fields["X"][:, k].ravel() for k in range(3))
    v = fields["U"][:, 1].ravel()
    psi, c = fields["S2"].ravel(), fields["S3"].ravel()

    w = 1.0 - psi
    xc, yc, zc = ((w*q).sum()/w.sum() for q in (x, y, z))
    r = np.hypot(x - xc, z - zc)
    u_in = -v[y > y.max() - 1.0].mean()

    def profile(coord, dist, lo, hi, width, *values):
        """Bin centres and the means of `values` over the points nearest (in
        `dist`) to the axis or plane in each bin of `coord` in [lo, hi)."""
        near = (dist < (args.radius or 0.2)) & (coord >= lo) & (coord < hi)
        k = np.floor((coord[near] - lo)/width).astype(int)
        dn, vals = dist[near], [q[near] for q in values]
        rows, widest = [], 0.0
        for b in np.unique(k):
            m = k == b
            # The mesh changes along the axis, so take the nearest points of each bin.
            m &= dn <= (args.radius or 1.5*dn[m].min() + 1e-9)
            widest = max(widest, dn[m].max())
            rows.append((lo + (b + 0.5)*width, *(q[m].mean() for q in vals)))
        return rows, widest

    def crossing(rows):
        """First place where psi (the second column) crosses 0.5, linearly
        interpolated between bins."""
        for (a, pa), (b, pb) in zip(rows, rows[1:]):
            if (pa - 0.5)*(pb - 0.5) <= 0 and pa != pb:
                return a + (0.5 - pa)*(b - a)/(pb - pa)
        raise SystemExit("psi does not cross 0.5: is there a bubble near the centroid?")

    axis, _ = profile(y, r, yc - 1.5, yc + 1.5, 0.01, psi)
    rear = crossing(axis)
    front = crossing(axis[::-1])
    equator, _ = profile(r, np.abs(y - yc), 0.0, 1.5, 0.01, psi)
    radius = crossing(equator)
    print(f"t = {time:.4f}: rear {rear:.4f}, front {front:.4f}, equatorial radius {radius:.4f}, "
          f"aspect ratio {radius/(0.5*(front - rear)):.2f}; inflow speed {u_in:.4f}")

    rows, widest = profile(y, r, rear - 3.0, rear + 0.05, 0.05, v, psi, c)
    nbins = int(round(3.05/0.05))
    if len(rows) < 0.5*nbins:
        raise SystemExit(f"only {len(rows)} of {nbins} bins have points near the axis: increase --radius")
    wake = [(rear - yy, vv/u_in, pp, cc) for yy, vv, pp, cc in rows[::-1]]
    print(f"  axis profile (points up to {widest:.4f} from the axis)")
    print("  below rear   v/u_in     c")
    for want in (0.05, 0.1, 0.2, 0.3, 0.5, 0.75, 1.0, 1.5, 2.0, 2.5):
        d, vv, pp, cc = min(wake, key=lambda row: abs(row[0] - want))
        print(f"  {d:9.3f}  {vv:+8.4f}  {cc:6.3f}")
    liquid = [(d, vv) for d, vv, pp, cc in wake if pp > 0.9]
    reverse = [d for d, vv in liquid if vv > 0]
    if reverse:
        print(f"  standing eddy: liquid moves up towards the bubble up to {max(reverse):.3f} D below it")
    else:
        print("  no standing eddy: the liquid on the axis moves away from the bubble everywhere below it")
    for frac in (0.02, 0.05, 0.1):
        slow = [d for d, vv in liquid if abs(vv) < frac]
        if slow:
            print(f"  |v| < {frac:g} u_in up to {max(slow):.3f} D below the rear")


if __name__ == "__main__":
    main()
