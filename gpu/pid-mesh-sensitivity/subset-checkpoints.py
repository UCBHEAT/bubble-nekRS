#!/usr/bin/env python3
"""Cut a run's checkpoints down to the elements animate.py draws, for meshes
too large for ParaView's Nek5000 reader.

Usage: subset-checkpoints.py [--margin 0.2] [--jobs 8] <run dir> <out dir>

ParaView 5.13's Nek5000 reader sizes its arrays with 32-bit integers and
aborts (std::bad_array_new_length) once 3 x elements x points per element
exceeds 2^31, above 1,398,101 elements at N = 7: on 0.132x (2,696,987
elements) of the Sc = 100 study, not yet on 0.198x (1,370,911). animate.py
needs only the elements that the z = 0 slice of its frames crosses and the
elements holding the bubble, apart from two of its bubble statistics, the
liquid mean velocity (u_l, v_l, w_l) and the largest liquid speed, which
span the whole domain.

For every checkpoint of <run dir>, read as animate.py reads them (each job's
part<N>/ with its own element order and the coordinates in its first file),
this writes a checkpoint with only the elements whose bounding box meets the
plane z = 0 or the box around the gas (psi < 0.99) of the part's first file,
grown by --margin bubble diameters, in the same part<N>/ layout under
<out dir>. From the full checkpoint it also computes those two whole-domain
statistics, with the point weights of animate.py's post, into
<out dir>/full-domain.csv, which animate.py merge puts into post.csv; and it
checks that the gas (1 - psi) a subset leaves out is negligible, so that the
other bubble statistics of the subset are those of the run: the margin must
span the tail of the psi profile, several interface widths eps, and the
bubble's changes of shape over a part (the Bo = 4 pilot's Nk27, eps = 0.032
and a wobbling bubble, needed 0.4; then its post.csv matched that of the
whole mesh to its 6 digits, but for 1e-7 in the centroid and gas means). <out dir> gets
links to the run's .par and data*.csv files, so that animate.py and
submit-animate-frontier.sh run on it as on the run (put it next to the run,
whose neighbours set the chart ranges), and the frames show what they would
show of the whole mesh.

Run it on a compute node: it reads every checkpoint whole (66 GB each on
0.132x), a part per process.
"""
import argparse
import csv
import glob
import multiprocessing
import os
import re
import sys

import numpy as np

CHUNK = 50000  # elements per read
HEADER = 132


def header(path):
    """The header of a Nek5000 #std FP32 file, with the offset and number of
    components of each field and of its per-element bounds (min and max of
    each component, which Nek5000 writes after the fields)."""
    with open(path, "rb") as f:
        raw = f.read(HEADER).decode()
        marker = np.frombuffer(f.read(4), dtype="<f4")[0]
    if abs(marker - 6.54321) > 1e-5:
        raise ValueError(f"{path}: not little-endian FP32 #std")
    words = raw.split()
    nxyz = int(words[2])*int(words[3])*int(words[4])
    nel, nelgt, code = int(words[5]), int(words[6]), words[11]
    if nel != nelgt:
        raise ValueError(f"{path}: one file per checkpoint expected")
    names, i = [], 0
    while i < len(code):
        if code[i] in "XU":
            names.append((code[i], 3))
            i += 1
        elif code[i] in "PT":
            names.append((code[i], 1))
            i += 1
        elif code[i] == "S":
            names += [(f"S{k + 1:02d}", 1) for k in range(int(code[i + 1:i + 3]))]
            i += 3
        else:
            raise ValueError(f"{path}: unknown field {code[i]} in {code}")
    blocks, offset = [], HEADER + 4 + 4*nel
    for name, ncomp in names:
        blocks.append([name, offset, ncomp, None])
        offset += 4*ncomp*nel*nxyz
    size = os.path.getsize(path)
    if size == offset + sum(8*ncomp*nel for _, ncomp in names):
        for block in blocks:
            block[3] = offset
            offset += 8*block[2]*nel
    if offset != size:
        raise ValueError(f"{path}: size {size}, header says {offset}")
    return {"raw": raw, "nxyz": nxyz, "n": int(words[2]), "nel": nel,
            "time": float(words[7]), "blocks": blocks}


def field(path, h, name, bounds=False):
    """A field of the file as a read-only array, (nel, 3, nxyz) for vectors;
    with bounds, its per-element bounds, (nel, 2 x components)."""
    for block, offset, ncomp, bounds_offset in h["blocks"]:
        if block == name:
            if bounds:
                return np.memmap(path, dtype="<f4", mode="r", offset=bounds_offset, shape=(h["nel"], 2*ncomp))
            shape = (h["nel"], 3, h["nxyz"]) if ncomp == 3 else (h["nel"], h["nxyz"])
            return np.memmap(path, dtype="<f4", mode="r", offset=offset, shape=shape)
    raise ValueError(f"{path}: no field {name}")


def weights_1d(points):
    """Trapezoid weights of the points along an element edge, (elements, 3,
    points): half the sum of the adjacent subcell widths, as animate.py's
    point_weights spreads each linear subcell's volume over its corners."""
    d = np.linalg.norm(np.diff(points, axis=2), axis=1)
    w = np.zeros(d.shape[:1] + (d.shape[1] + 1,))
    w[:, :-1] += 0.5*d
    w[:, 1:] += 0.5*d
    return w


def select(first, h, margin):
    """Elements meeting the plane z = 0 or the gas box of the part's first file
    grown by margin; the gas box; the edge weights of every element along its
    three point index directions (which gmsh2nek's elements do not align with
    x, y and z in any fixed order; the elements are boxes, as point_weights
    assumes)."""
    n = h["n"]
    X = field(first, h, "X")
    psi = field(first, h, "S02")
    lo, hi = np.empty((h["nel"], 3)), np.empty((h["nel"], 3))
    w1 = np.empty((h["nel"], 3, n))
    gas_lo, gas_hi = np.full(3, np.inf), np.full(3, -np.inf)
    for a in range(0, h["nel"], CHUNK):
        b = min(a + CHUNK, h["nel"])
        x = np.asarray(X[a:b], dtype=np.float64)
        lo[a:b], hi[a:b] = x.min(axis=2), x.max(axis=2)
        p = x.reshape(b - a, 3, n, n, n)  # [element, component, k, j, i], i fastest
        w1[a:b, 0] = weights_1d(p[:, :, 0, 0, :])
        w1[a:b, 1] = weights_1d(p[:, :, 0, :, 0])
        w1[a:b, 2] = weights_1d(p[:, :, :, 0, 0])
        gas = np.asarray(psi[a:b]) < 0.99
        if gas.any():
            pts = x.transpose(0, 2, 1)[gas]
            gas_lo, gas_hi = np.minimum(gas_lo, pts.min(axis=0)), np.maximum(gas_hi, pts.max(axis=0))
    box_lo, box_hi = gas_lo - margin, gas_hi + margin
    in_box = ((hi >= box_lo) & (lo <= box_hi)).all(axis=1)
    on_plane = (lo[:, 2] <= 0.0) & (hi[:, 2] >= 0.0)
    return in_box | on_plane, (gas_lo, gas_hi), w1


def process(path, h, keep, w1, out):
    """Write the file with only the elements keep (numbered 1.. in their
    order), reading it once, and return the whole-domain liquid mean velocity
    and largest liquid speed (as animate.py's post computes them) and the
    gas outside keep."""
    nxyz, nel = h["nxyz"], h["nel"]
    m = int(keep.sum())
    before = np.concatenate([[0], np.cumsum(keep)])  # kept elements before each element
    # The element counts are the two i10 fields after "#std w nn nn nn".
    raw = h["raw"]
    match = re.match(r"(#std\s+\d+\s+\d+\s+\d+\s+\d+)(\s+\d+)(\s+\d+)", raw)
    head = match.group(1) + f"{m:>{len(match.group(2))}d}" + f"{m:>{len(match.group(3))}d}" + raw[match.end():]
    assert len(head) == HEADER
    out_offset, offset = {}, HEADER + 4 + 4*m
    for name, _, ncomp, _ in h["blocks"]:
        out_offset[name] = offset
        offset += 4*ncomp*m*nxyz
    data = {name: field(path, h, name) for name, _, _, _ in h["blocks"]}

    sums = np.zeros(3)
    liquid = gas = gas_out = speed_max = 0.0
    tmp = out + ".tmp"
    with open(tmp, "wb") as f:
        f.write(head.encode())
        f.write(np.array([6.54321], dtype="<f4").tobytes())
        f.write(np.arange(1, m + 1, dtype="<i4").tobytes())
        for a in range(0, nel, CHUNK):
            b = min(a + CHUNK, nel)
            k = keep[a:b]
            chunk = {}
            for name, _, ncomp, _ in h["blocks"]:
                if name in ("U", "S02") or k.any():
                    chunk[name] = np.asarray(data[name][a:b])
                    if k.any():
                        f.seek(out_offset[name] + 4*ncomp*nxyz*int(before[a]))
                        f.write(np.ascontiguousarray(chunk[name][k]).tobytes())
            w = w1[a:b]
            W = (w[:, 2, :, None, None]*w[:, 1, None, :, None]*w[:, 0, None, None, :]).reshape(b - a, nxyz)
            psi = np.clip(chunk["S02"].astype(np.float64), 0.0, 1.0)
            u = chunk["U"].astype(np.float64)
            wl = W*psi
            liquid += wl.sum()
            sums += np.einsum("ep,ecp->c", wl, u)
            wg = W*(1.0 - psi)
            gas += wg.sum()
            gas_out += wg[~k].sum()
            liq = psi > 0.95
            if liq.any():
                speed_max = max(speed_max, np.sqrt((u*u).sum(axis=1))[liq].max())
        # The kept elements' bounds, in the same order as the fields.
        f.seek(offset)
        for name, _, _, bounds_offset in h["blocks"]:
            if bounds_offset is not None:
                f.write(np.ascontiguousarray(np.asarray(field(path, h, name, bounds=True))[keep]).tobytes())
    os.rename(tmp, out)
    ul = sums/liquid
    return m, {"Time": h["time"], "u_l": ul[0], "v_l": ul[1], "w_l": ul[2],
               "max_speed_liquid": speed_max, "gas_volume": gas, "gas_outside_subset": gas_out}


def do_part(args):
    part, out_part, margin = args
    files = sorted(glob.glob(os.path.join(part, "*0.f[0-9][0-9][0-9][0-9][0-9]")))
    if not files:
        return []
    h0 = header(files[0])
    keep, (gas_lo, gas_hi), w1 = select(files[0], h0, margin)
    os.makedirs(out_part, exist_ok=True)
    rows = []
    for path in files:
        h = header(path)
        if h["nel"] != h0["nel"]:
            raise ValueError(f"{path}: element count differs from {files[0]}")
        m, row = process(path, h, keep, w1, os.path.join(out_part, os.path.basename(path)))
        rows.append(row)
        print(f"{path}: t = {h['time']:.4f}, kept {m} of {h['nel']} elements, gas box "
              f"{np.round(gas_lo, 3)} to {np.round(gas_hi, 3)}, gas {row['gas_volume']:.6f} "
              f"({row['gas_outside_subset']:.2e} outside), u_l {row['u_l']:.6f} {row['v_l']:.6f} "
              f"{row['w_l']:.6f}, liquid speed max {row['max_speed_liquid']:.4f}", flush=True)
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--margin", type=float, default=0.2,
                        help="growth of the gas box, in bubble diameters (default 0.2)")
    parser.add_argument("--jobs", type=int, default=8, help="parts read at once (default 8)")
    parser.add_argument("run")
    parser.add_argument("out")
    args = parser.parse_args()

    run, out = os.path.abspath(args.run), os.path.abspath(args.out)
    parts = sorted(glob.glob(os.path.join(run, "part[0-9]*")), key=lambda d: int(os.path.basename(d)[4:]))
    targets = [(p, os.path.join(out, os.path.basename(p)), args.margin) for p in parts] or \
              [(run, out, args.margin)]
    os.makedirs(out, exist_ok=True)
    for path in glob.glob(os.path.join(run, "*.par")) + glob.glob(os.path.join(run, "data*.csv")):
        link = os.path.join(out, os.path.basename(path))
        if not os.path.lexists(link):
            os.symlink(path, link)

    # The workers raise ordinary exceptions: a SystemExit in a pool worker
    # would hang pool.map.
    try:
        with multiprocessing.Pool(min(args.jobs, len(targets))) as pool:
            rows = [row for part_rows in pool.map(do_part, targets) for row in part_rows]
    except ValueError as e:
        sys.exit(str(e))
    rows.sort(key=lambda r: r["Time"])
    worst = max(r["gas_outside_subset"]/r["gas_volume"] for r in rows)
    print(f"largest fraction of the gas outside a subset: {worst:.2e}")
    if worst > 1e-6:
        sys.exit("the subsets miss gas: grow --margin")
    with open(os.path.join(out, "full-domain.csv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=["Time", "u_l", "v_l", "w_l", "max_speed_liquid"],
                           extrasaction="ignore", lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: repr(float(v)) for k, v in r.items()})


if __name__ == "__main__":
    main()
