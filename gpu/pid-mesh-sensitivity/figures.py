#!/usr/bin/env python3
"""README figures of the meshes and of the concentration wake.

Usage: figures.py [--runs <dir with the run directories>] [--out doc]

doc/mesh.png: the z = 0 plane of the five meshes (element outlines, from
each run's bubble.plan.json) with the initial bubble (red circle), and a
zoom of the 1x mesh around the bubble at t = 30 with its GLL points and the
psi = 0.5 interface.

doc/wake.png: c on the z = 0 plane of the 0.444x run at t = 30, over the
whole domain and around the bubble, with the interface and the bubble-frame
streamlines.

The slices are interpolated to z = 0 with each element's own GLL basis, from
the last checkpoint of the run (with the coordinates of that job's first
checkpoint). The meshes are tensor products of axis-aligned elements, so
z = 0 crosses one layer of elements (the bubble centre lies inside an
element, see the README).
"""
import argparse
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PolyCollection
from matplotlib.tri import LinearTriInterpolator, Triangulation
from numpy.polynomial import legendre

RUNS = ["2.25x", "1.5x", "1x", "0.667x", "0.444x"]


def gll_nodes(n):
    """The n Gauss-Lobatto-Legendre nodes on [-1, 1] (N = n - 1)."""
    p = np.zeros(n)
    p[-1] = 1.0
    return np.concatenate(([-1.0], np.sort(legendre.legroots(legendre.legder(p))), [1.0]))


def lagrange(nodes, x):
    """Values at x of the Lagrange polynomials on the nodes."""
    out = np.ones(len(nodes))
    for k, xk in enumerate(nodes):
        for m, xm in enumerate(nodes):
            if m != k:
                out[k] *= (x - xm)/(xk - xm)
    return out


def header(path):
    with open(path, "rb") as f:
        words = f.read(132).decode().split()
    n = int(words[2])
    return n, int(words[5]), words[11], float(words[7])


def last_checkpoint(run):
    """The run's last checkpoint and the first checkpoint of its job (which
    holds the coordinates), from the time-ordered links of finish-runs.sh."""
    last = os.path.realpath(sorted(glob.glob(os.path.join(run, "bubble0.f[0-9]*")))[-1])
    return last, os.path.join(os.path.dirname(last), "bubble0.f00000")


def z_slice(ckpt, coords, z0=0.0):
    """x, y, triangles and the fields (u, v, psi, c) on the plane z = z0."""
    n, nel, code, time = header(ckpt)
    nx, _, code_x, _ = header(coords)
    assert code_x.startswith("X") and nx == n
    npt = n**3
    head = 136 + 4*nel
    X = np.memmap(coords, dtype="<f4", mode="r", offset=head, shape=(nel, 3, npt))
    zc = np.asarray(X[:, 2, :])
    hit = np.where((zc.min(axis=1) <= z0) & (zc.max(axis=1) > z0))[0]
    # Field blocks of the checkpoint, in nel*npt units: [X,] U (3), P, S01-S03.
    base = 3 if code.startswith("X") else 0
    block = 4*nel*npt

    def field(index, comp=None):
        if comp is None:
            a = np.memmap(ckpt, dtype="<f4", mode="r", offset=head + index*block, shape=(nel, npt))
            return np.asarray(a[hit])
        a = np.memmap(ckpt, dtype="<f4", mode="r", offset=head + index*block, shape=(nel, 3, npt))
        return np.asarray(a[hit, comp])

    arrays = {"x": np.asarray(X[hit, 0]), "y": np.asarray(X[hit, 1]), "z": zc[hit],
              "u": field(base, 0), "v": field(base, 1),
              "psi": field(base + 5), "c": field(base + 6)}
    nodes = gll_nodes(n)
    out = {k: np.empty((len(hit), n, n)) for k in arrays if k != "z"}
    for e in range(len(hit)):
        z = arrays["z"][e].reshape(n, n, n)          # [k][j][i]
        spans = [abs(z[0, 0, -1] - z[0, 0, 0]), abs(z[0, -1, 0] - z[0, 0, 0]), abs(z[-1, 0, 0] - z[0, 0, 0])]
        axis = 2 - int(np.argmax(spans))             # array axis along which z varies
        line = np.moveaxis(z, axis, 0)[:, 0, 0]
        xi = 2.0*(z0 - line[0])/(line[-1] - line[0]) - 1.0
        w = lagrange(nodes, xi)
        for k in out:
            out[k][e] = np.tensordot(w, np.moveaxis(arrays[k][e].reshape(n, n, n), axis, 0), axes=1)
    # Two triangles per GLL sub-quad, within each element.
    idx = np.arange(n*n).reshape(n, n)
    quad = np.stack([idx[:-1, :-1], idx[:-1, 1:], idx[1:, 1:], idx[1:, :-1]], axis=-1).reshape(-1, 4)
    one = np.concatenate([quad[:, [0, 1, 2]], quad[:, [0, 2, 3]]])
    tris = (one[None, :, :] + (n*n*np.arange(len(hit)))[:, None, None]).reshape(-1, 3)
    flat = {k: v.reshape(-1) for k, v in out.items()}
    # Neighbouring elements share their edge points; merge them (the fields
    # are continuous there), since the triangle finder needs unique points.
    _, inv = np.unique(np.round(np.column_stack([flat["x"], flat["y"]]), 6), axis=0, return_inverse=True)
    inv = inv.ravel()
    count = np.bincount(inv)
    flat = {k: np.bincount(inv, weights=v)/count for k, v in flat.items()}
    return Triangulation(flat["x"], flat["y"], inv[tris]), flat, time


def mesh_lines(run):
    """Element node coordinates along x (= z) and y, from bubble.plan.json."""
    with open(os.path.join(run, "bubble.plan.json")) as f:
        layout = json.load(f)["layout"]
    xs = -3.0 + np.concatenate(([0.0], np.cumsum(layout["x_and_z"]["sizes"])))
    ys = -10.0 + np.concatenate(([0.0], np.cumsum(layout["y"]["sizes"])))
    return xs, ys, layout


def mesh_figure(runs_dir, out):
    lines = {r: mesh_lines(os.path.join(runs_dir, r)) for r in RUNS}
    fig = plt.figure(figsize=(13.5, 6.0))
    left = fig.add_gridspec(1, 5, wspace=0.08, left=0.04, right=0.645, top=0.86, bottom=0.07)
    right = fig.add_gridspec(1, 1, left=0.695, right=0.99, top=0.86, bottom=0.07)
    for k, r in enumerate(RUNS):
        ax = fig.add_subplot(left[0, k])
        xs, ys, layout = lines[r]
        polys = [[(xs[i], ys[j]), (xs[i + 1], ys[j]), (xs[i + 1], ys[j + 1]), (xs[i], ys[j + 1])]
                 for i in range(len(xs) - 1) for j in range(len(ys) - 1)]
        ax.add_collection(PolyCollection(polys, facecolors="none", edgecolors="0.3", linewidths=0.3))
        ax.add_patch(plt.Circle((0, 0), 0.5, fill=False, color="#d62728", lw=1.0))
        ax.set_xlim(-3, 3)
        ax.set_ylim(-10, 5)
        ax.set_aspect("equal")
        ax.set_xticks([-2, 0, 2])
        ax.set_yticks([-10, -5, 0, 5])
        if k:
            ax.set_yticklabels([])
        n_xz, n_y = layout["x_and_z"]["n_elements"], layout["y"]["n_elements"]
        ax.set_title("{}\n{}x{}x{} = {}".format(r, n_xz, n_y, n_xz, n_xz*n_y*n_xz), fontsize=10)
        ax.tick_params(labelsize=9)
    fig.text(0.04, 0.95, "z = 0 plane of the five meshes: inflow (c = 1) at the top, outflow at "
             "the bottom, x and z periodic; red: initial bubble (D = 1)", fontsize=10, ha="left")

    # Zoom: the 1x mesh around the bubble at t = 30.
    run = os.path.join(runs_dir, "1x")
    tri, f, time = z_slice(*last_checkpoint(run))
    ax = fig.add_subplot(right[0, 0])
    xs, ys, layout = lines["1x"]
    xlim, ylim = (-1.15, 1.15), (-0.95, 0.85)
    nodes = gll_nodes(8)
    gx = np.concatenate([xs[i] + 0.5*(nodes + 1)*(xs[i + 1] - xs[i]) for i in range(len(xs) - 1)])
    gy = np.concatenate([ys[j] + 0.5*(nodes + 1)*(ys[j + 1] - ys[j]) for j in range(len(ys) - 1)])
    gx, gy = gx[(gx >= xlim[0]) & (gx <= xlim[1])], gy[(gy >= ylim[0]) & (gy <= ylim[1])]
    GX, GY = np.meshgrid(gx, gy)
    ax.plot(GX.ravel(), GY.ravel(), ".", color="0.55", ms=1.6, zorder=1)
    for xv in xs[(xs >= xlim[0]) & (xs <= xlim[1])]:
        ax.axvline(xv, color="0.15", lw=0.8, zorder=2)
    for yv in ys[(ys >= ylim[0]) & (ys <= ylim[1])]:
        ax.axhline(yv, color="0.15", lw=0.8, zorder=2)
    ax.tricontour(tri, f["psi"], levels=[0.5], colors="#ff7f0e", linewidths=1.8, zorder=3)
    ax.add_patch(plt.Circle((0, 0), 0.5, fill=False, color="#d62728", lw=1.0, ls="--", zorder=3))
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_aspect("equal")
    ax.tick_params(labelsize=9)
    ax.set_title("1x around the bubble at t = {:.0f}: elements, GLL points\n(N = 7), "
                 "interface psi = 0.5 (orange), initial bubble (red)".format(time), fontsize=10)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(out)


def wake_figure(runs_dir, out):
    run = os.path.join(runs_dir, "0.444x")
    tri, f, time = z_slice(*last_checkpoint(run))
    gas = f["psi"] < 0.5
    cmap = plt.get_cmap("Blues_r")
    levels = np.linspace(0.0, 1.0, 21)
    levels[-1] = 1.001

    fig = plt.figure(figsize=(10.5, 7.6))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.55, 0.06], wspace=0.12,
                          left=0.06, right=0.93, top=0.9, bottom=0.07)
    for k, (xlim, ylim) in enumerate([((-3, 3), (-10, 5)), ((-1.3, 1.3), (-2.8, 1.0))]):
        ax = fig.add_subplot(gs[0, k])
        cs = ax.tricontourf(tri, np.clip(f["c"], 0.0, 1.0), levels=levels, cmap=cmap)
        ax.tricontourf(tri, f["psi"], levels=[-1.0, 0.5], colors=["0.85"])
        ax.tricontour(tri, f["psi"], levels=[0.5], colors="k", linewidths=1.0)
        if k == 1:
            # Bubble-frame streamlines (the liquid enters at the top).
            gx = np.linspace(*xlim, 220)
            gy = np.linspace(*ylim, 320)
            GX, GY = np.meshgrid(gx, gy)
            U = LinearTriInterpolator(tri, f["u"])(GX, GY)
            V = LinearTriInterpolator(tri, f["v"])(GX, GY)
            inside = LinearTriInterpolator(tri, f["psi"])(GX, GY) < 0.5
            U, V = np.ma.masked_where(inside, U), np.ma.masked_where(inside, V)
            ax.streamplot(gx, gy, U, V, density=1.4, color="w", linewidth=0.6, arrowsize=0.7)
        ax.set_xlim(*xlim)
        ax.set_ylim(*ylim)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=9)
        ax.set_xlabel("x / D", fontsize=10)
        if k == 0:
            ax.set_ylabel("y / D", fontsize=10)
            ax.add_patch(plt.Rectangle((-1.3, -2.8), 2.6, 3.8, fill=False, ec="k", lw=0.8, ls="--"))
    cax = fig.add_subplot(gs[0, 2])
    cb = fig.colorbar(cs, cax=cax, ticks=np.linspace(0, 1, 6))
    cb.set_label("c (1 = saturated inflow, 0 in the gas)", fontsize=10)
    cb.ax.tick_params(labelsize=9)
    fig.suptitle("0.444x mesh, z = 0 plane, t = {:.0f}: the concentration wake (gray: gas; right: "
                 "zoom with bubble-frame streamlines)".format(time), fontsize=11)
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(out)


def main():
    here = os.path.dirname(os.path.abspath(__file__))
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--runs", default="/lus/eagle/projects/nek-vf/benl/answinter26/Sc1-pid",
                        help="directory holding the run directories")
    parser.add_argument("--out", default=os.path.join(here, "doc"))
    args = parser.parse_args()
    os.makedirs(args.out, exist_ok=True)
    mesh_figure(args.runs, os.path.join(args.out, "mesh.png"))
    wake_figure(args.runs, os.path.join(args.out, "wake.png"))


if __name__ == "__main__":
    main()
