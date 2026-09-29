r"""Hex mesh for the 3D PID-centred rising bubble (nekRS-LS), built from gmsh transfinite blocks.

The bubble (D = 1) sits at the origin, y is vertical, and in the PID frame the
liquid enters at the top (inflow, y = +Lu) and leaves at the bottom (outflow,
y = -Ld), so the wake is below the bubble. x and z are periodic with period W.

Layout (z = 0 cut; x and z are treated identically):

      y = +Lu  +------+-------------+------+   inflow  (boundary ID 1)
               |      :  upstream   :      |   sizes grow to <= h_far
               |      :             :      |
         +b    +------+-------------+------+
               | side |  core block | side |   |x|,|y|,|z| <= b around the bubble (o),
               |      |     (o)     |      |   see below
         -b    +------+-------------+------+
               |      | wake column |      |   cross-section = the core block face's;
               |      | h_wake_near |      |   streamwise: h_wake_near, linear ramp
               |      |   ramp      |      |   to h_wake_far, h_wake_far, then grow
               |      | h_wake_far  |      |   to h_far
               |      :             :      |
      y = -Ld  +------+-------------+------+   outflow (boundary ID 2)
             -W/2    -b             +b    +W/2   xmin/xmax, zmin/zmax periodic

All blocks are axis-aligned and conforming, so outside the core block the mesh
is a global tensor product of three 1-D node distributions (x = z, and y): the
core block's grid lines run through every block to the boundaries. Along each
axis the elements are uniform across the core block and grow away from it by
at most a factor `growth` per element (checked: axis_ratio_max), never
exceeding the local size target (the wake schedule below the core block,
h_wake_cross beside it out to rho = r_wake, h_far elsewhere). Each
distribution is split into runs of constant size ratio (uniform or
geometric); every run becomes one segment of the block grid, meshed with a
gmsh transfinite "Progression" curve, so the transfinite blocks reproduce the
designed distribution exactly (checked after meshing). The x and z
distributions are mirror-symmetric, which makes the periodic faces match node
for node.

core = "cartesian": the core block is a uniform fine cube of half-width
b = fine_half with cells of h = min(h_iface, h_wake_cross, h_wake_near); by
default fine_half is the smallest multiple of h that covers the interface
shell r <= r_iface_out + iface_margin. Every element is an axis-aligned brick.

core = "ogrid": the core block [-a, a]^3 (b = a) has n x n cells of
h_core = 2a/n = frac h_wake_cross <= h_wake_cross on each face (frac from
BOX_FRACTIONS), so its faces continue the wake column and the side blocks,
and it is filled with an O-grid (one octant, y = 0 cut):

      +---------------+  core box face, flat, n cells of h_core = 2a/n
      |\  layer K     |  levels: closed "bulged cubes" (6 spherical caps
      | \  ...        |  through 8 corners); level 1 lies just outside the
      |  +--------+   |  interface shell (face-centre radius r_iface_out +
      |  | layer 1|   |  iface_margin), the last level is the core box itself
      |  | +----+ |   |
      |  | |cube| |   |  flat inner cube, n^3 cells of 2c/n, holding the
      |  | |(o) | |   |  bubble; cube cells and layer 1 are ~isotropic and
      +--+-+----+-+---+  <= h_iface

The O-grid decouples the interface cell size from h_core (the wake-column and
side spacing), which a tensor product cannot do; it saves elements for
h_iface below about 0.16 (at 0.2 the cartesian core is cheaper). The layout is
derived from the targets without trial meshing: the level nodes are computed
exactly (level_nodes, the construction gmsh uses), so the edges of any
candidate layer are known before meshing. ogrid_candidate derives the layout
for a cell count n and box fraction frac: the inner cube and the flattest
admissible first layer meet the interface targets, the outer layers grow by
<= growth, at the face centres and along the corner rays, under the largest
thickness that keeps the tilted radial edges below the bubble within the wake
targets, and the side and wake cells outside the box start from the last
layer's face-centre thickness; a candidate with an edge below
edge_min_fraction h_iface or an axis size ratio above growth is ruled out.
ogrid_layout keeps the candidate with the fewest elements in the whole mesh.

After meshing, the element types and count and the designed node positions
are checked (check_nodes), and every zone target (check_mesh.py); a failed
target raises an error naming the parameter to change. Output: NAME.msh (msh
2.2 binary, Hexahedron 27) with the physical groups, in this order (gmsh2nek
numbers boundaries 1..6 in order): inflow, outflow, xmin, xmax, zmin, zmax,
and the volume group fluid; and NAME.plan.json (parameters, layout, element
count, checks and the derived nekRS-LS settings).

usage: python generate_mesh.py [--preset NAME] [--KEY VALUE ...]
       e.g. --preset fine, --preset ogrid --W 8, --h_iface 0.15
"""
import argparse
import itertools
import json
import math
import time

import gmsh
import numpy as np

import check_mesh

# Parameters (lengths in bubble diameters, D = 1). Presets below and
# command-line options (--KEY VALUE) override them.
PARAMS = {
    # Domain: x, z in [-W/2, W/2] (periodic), y in [-Ld, +Lu].
    "W": 6.0,
    "Lu": 5.0,
    "Ld": 10.0,
    # Core block around the bubble: "cartesian" (uniform fine cube) or
    # "ogrid" (inner cube and O-grid layers in a coarser box).
    "core": "cartesian",
    # Interface zone: every element overlapping the shell
    # r_iface_in <= r <= r_iface_out has edges <= h_iface and an edge aspect
    # ratio <= aspect_iface. The fine cube (cartesian) or the first O-grid
    # level (ogrid) lies at least iface_margin outside r_iface_out.
    "r_iface_in": 0.3,
    "r_iface_out": 0.75,
    "h_iface": 0.2,
    "aspect_iface": 2.0,
    "iface_margin": 0.01,
    # cartesian only: half-width of the fine cube (0 = derived from the
    # interface zone; an explicit value must cover r_iface_out + iface_margin).
    "fine_half": 0.0,
    # Wake zone: every other element overlapping the column rho <= r_wake,
    # y_wake_end <= y <= y_wake_top has cross-stream edges <= h_wake_cross
    # and streamwise edges <= h_wake_near down to y_wake_near, then a linear
    # ramp to h_wake_far at y_wake_ramp, then h_wake_far down to y_wake_end.
    # In ogrid mode the core box face cells are h_core = frac h_wake_cross
    # (frac from BOX_FRACTIONS).
    "r_wake": 0.75,
    "y_wake_top": -0.5,
    "h_wake_cross": 0.3,
    "h_wake_near": 0.25,
    "y_wake_near": -1.5,
    "h_wake_far": 0.5,
    "y_wake_ramp": -5.0,
    "y_wake_end": -6.0,
    # Far field: largest element edge, and the largest size ratio of
    # neighbouring elements along an axis (also between O-grid layers; at
    # most thickness_ratio_max).
    "h_far": 1.25,
    "growth": 1.5,
    # Further targets checked after meshing (see check_mesh.py): neighbour
    # h_J ratio, face-normal thickness ratio, gmsh minSJ, smallest edge as a
    # fraction of h_iface.
    "hJ_ratio_max": 1.5,
    "thickness_ratio_max": 2.0,
    "minSJ_min": 0.5,
    "edge_min_fraction": 0.5,
    # nekRS-LS settings in NAME.plan.json are derived for N = 7 and N_target
    # (Weber number and liquid/gas density ratio of the case).
    "N_target": 7,
    "We": 2.667,
    "rho_ratio": 40.0,
    # Output mesh name (NAME.msh, NAME.plan.json) and gmsh threads.
    "name": "bubble",
    "threads": 4,
    # Safety limits: refuse a layout with more elements than this (checked
    # before any large array is built), and cap the process's address space
    # (GiB, 0: none; see check_mesh.limit_memory), so a bad parameter set fails
    # with an error instead of exhausting the machine's memory.
    "max_elements": 200000,
    "max_memory_gb": 8.0,
}
CORES = ("cartesian", "ogrid")
# Parameters that must be positive (all float parameters must be finite).
POSITIVE = ("W", "Lu", "Ld", "r_iface_out", "h_iface", "aspect_iface", "iface_margin", "r_wake", "h_wake_cross",
            "h_wake_near", "h_wake_far", "h_far", "hJ_ratio_max", "thickness_ratio_max", "We", "rho_ratio",
            "max_elements")

# Named parameter sets (--preset NAME); command-line options override them.
PRESETS = {
    "default": {},                                          # 0.2 cubes around the bubble
    "fine": {"h_iface": 0.125, "h_wake_cross": 0.2},        # 0.125 cubes (also in the wake column)
    "large": {"W": 8.0, "Lu": 6.0, "Ld": 15.0},             # default resolution, larger domain
    "ogrid": {"core": "ogrid", "h_iface": 0.15, "h_wake_cross": 0.3},
    "ogrid-fine": {"core": "ogrid", "h_iface": 0.1, "h_wake_cross": 0.2},
}

# Graded sizes and size ratios are designed this much (relative) below their
# targets and growth: gmsh places the nodes of geometric Progression curves by
# numerical integration, accurate to ~1e-8, which would otherwise push sizes
# designed exactly at a target (e.g. h_wake_far), or ratios at growth =
# thickness_ratio_max, just above it.
SIZE_MARGIN = 1 - 1e-6
# Number of intervals of the table (WAKE_TABLE + 1 radii) of the largest
# wake-admissible O-grid layer thickness.
WAKE_TABLE = 16
# O-grid layout search (ogrid_layout): how many cell counts n to try (not
# counting those whose largest core box is too small), the core-box
# half-widths tried as fractions of n h_wake_cross / 2 (the fraction sets the
# wake-column cross-stream spacing), and the most layers accepted.
N_OGRID_TRIES = 3
BOX_FRACTIONS = (1.0, 0.95, 0.9, 0.85, 0.8, 0.75, 0.7)
MAX_OGRID_LAYERS = 30


def graded_sizes(length, h_prev, target, growth, max_n):
    """Fewest element sizes filling `length` outward from a face of the core block.

    `h_prev` is the size of the element on the other side of the face and
    `target(s)` the largest allowed size of an element starting at distance s
    from the face (non-decreasing in s). Each size is the largest one allowed
    by the target and by `growth` times the previous size (both times
    SIZE_MARGIN); this greedy march needs the fewest elements. Its n elements
    are then fitted exactly to `length` by capping all sizes at a common value
    c (found by bisection), which shrinks only the largest, outermost
    elements. Neighbour ratios stay <= growth, except next to the face when
    the target or c is below h_prev / growth. Raises ValueError if more than
    max_n elements are needed.
    """
    def march(cap, n=None):
        sizes, s, h = [], 0.0, h_prev
        while (s < length * (1 - 1e-12)) if n is None else (len(sizes) < n):
            h = min(SIZE_MARGIN * growth * h, SIZE_MARGIN * target(s), cap)
            sizes.append(h)
            s += h
            if len(sizes) > max_n:
                raise ValueError(f"more than {max(max_n, 0)} graded elements needed (the layout would exceed "
                                 "max_elements)")
        return np.array(sizes)

    free = march(math.inf)
    n = len(free)
    lo, hi = 0.0, free.max()
    for _ in range(200):
        c = 0.5 * (lo + hi)
        if march(c, n).sum() > length:
            hi = c
        else:
            lo = c
    sizes = march(hi, n)
    return sizes * (length / sizes.sum())


def wake_target(p, y0):
    """Allowed streamwise size of an element whose upper face lies at distance
    s below y = y0 (the bottom of the core block): the wake schedule, then h_far."""
    def target(s):
        y = y0 - s
        if y >= p["y_wake_near"]:
            return p["h_wake_near"]
        if y >= p["y_wake_ramp"]:
            t = (p["y_wake_near"] - y) / (p["y_wake_near"] - p["y_wake_ramp"])
            return p["h_wake_near"] + t * (p["h_wake_far"] - p["h_wake_near"])
        if y >= p["y_wake_end"]:
            return p["h_wake_far"]
        return p["h_far"]
    return target


def runs(sizes, rtol=1e-9):
    """Lengths of the fewest consecutive runs of constant size ratio (uniform
    or geometric) that make up `sizes`; each run can be meshed by one
    transfinite Progression curve. Ties are broken toward fewer elements in
    non-uniform runs, e.g. [0.5, 0.5], [0.75], [1, 1, 1] rather than
    [0.5, 0.5], [0.75, 1], [1, 1] (uniform curves are meshed exactly)."""
    def ratio_is(i, j, r):
        q = sizes[i + 1:j] / sizes[i:j - 1]
        return bool(np.all(np.abs(q - r) <= rtol * r))

    n = len(sizes)
    best = [(0, 0, [])] + [None] * n  # (runs, elements in non-uniform runs, run lengths)
    for j in range(1, n + 1):
        for i in range(j):
            if best[i] is None or not ratio_is(i, j, sizes[i + 1] / sizes[i] if j - i > 1 else 1.0):
                continue
            cand = (best[i][0] + 1, best[i][1] + (0 if ratio_is(i, j, 1.0) else j - i), best[i][2] + [j - i])
            if best[j] is None or cand[:2] < best[j][:2]:
                best[j] = cand
    return best[n][2]


class Axis:
    """One axis of the tensor-product mesh: the element sizes along [lo, hi]
    (n_core uniform elements across the core block |t| <= half, graded
    sides), the node coordinates, and the block breakpoints (ends of the
    constant-ratio runs); segment `core` is the core block.

    h_prev is the thickness of the core block's outermost element at its face
    centres (the first side element is at most growth * h_prev), target_lo
    and target_hi give the allowed element size as a function of the distance
    from the core block on each side (see graded_sizes), and the axis may
    have at most max_n elements (else ValueError).
    """

    def __init__(self, lo, hi, half, n_core, h_prev, target_lo, target_hi, growth, max_n):
        if not lo < -half < half < hi:
            raise ValueError(f"core block [-{half:g}, {half:g}] must lie inside [{lo:g}, {hi:g}]: "
                             "enlarge the domain (W, Lu, Ld)")
        below = graded_sizes(-half - lo, h_prev, target_lo, growth, max_n - n_core)  # outward from the core
        above = graded_sizes(hi - half, h_prev, target_hi, growth, max_n - n_core - len(below))
        self.sizes = np.concatenate([below[::-1], np.full(n_core, 2 * half / n_core), above])
        self.nodes = lo + np.concatenate([[0.0], np.cumsum(self.sizes)])
        self.nodes[-1] = hi
        r_below = runs(below)[::-1]
        self.cuts = np.cumsum([0] + r_below + [n_core] + runs(above))
        self.breaks = self.nodes[self.cuts]
        self.core = len(r_below)

    def segment(self, a):
        """(number of nodes, Progression coefficient) of block segment a."""
        i0, i1 = self.cuts[a], self.cuts[a + 1]
        coef = self.sizes[i0 + 1] / self.sizes[i0] if i1 - i0 > 1 else 1.0
        if abs(coef - 1) < 1e-9:
            coef = 1.0  # uniform: gmsh then places the nodes exactly
        return int(i1 - i0 + 1), float(coef)

    def ratio(self):
        """Largest size ratio of neighbouring elements."""
        q = self.sizes[1:] / self.sizes[:-1]
        return float(np.maximum(q, 1 / q).max())

    def plan(self):
        """JSON-ready description of the axis."""
        return {"n_elements": len(self.sizes), "sizes": self.sizes, "breakpoints": self.breaks,
                "segments": [self.segment(a) for a in range(len(self.breaks) - 1)], "size_ratio_max": self.ratio()}


def validate(p):
    """Reject inconsistent parameters before meshing."""
    if p["core"] not in CORES:
        raise ValueError(f"core must be one of {CORES}")
    for k, v in p.items():
        if isinstance(v, float) and not math.isfinite(v):
            raise ValueError(f"{k} = {v} must be finite")
    bad = [k for k in POSITIVE if not p[k] > 0]
    if bad:
        raise ValueError(f"{', '.join(bad)} must be positive")
    if not (p["N_target"] >= 2 and p["threads"] >= 0 and p["fine_half"] >= 0):
        raise ValueError("need N_target >= 2, threads >= 0 and fine_half >= 0")
    if not (p["h_wake_near"] <= p["h_wake_far"] <= p["h_far"] and p["h_wake_cross"] <= p["h_far"]
            and p["y_wake_top"] >= p["y_wake_near"] >= p["y_wake_ramp"] >= p["y_wake_end"] > -p["Ld"]):
        raise ValueError("wake schedule must coarsen downstream: h_wake_near <= h_wake_far <= h_far, "
                         "h_wake_cross <= h_far, y_wake_top >= y_wake_near >= y_wake_ramp >= y_wake_end > -Ld")
    if not (0 <= p["r_iface_in"] < p["r_iface_out"] and 1 < p["growth"] <= p["thickness_ratio_max"]):
        raise ValueError("need 0 <= r_iface_in < r_iface_out and 1 < growth <= thickness_ratio_max")
    if p["fine_half"] and p["core"] != "cartesian":
        raise ValueError("fine_half applies to core = cartesian only")


def targets(p):
    """The check_mesh targets of the parameters p."""
    return {k: p[k] for k in check_mesh.TARGETS}


def cartesian_layout(p):
    """Fine cube of the cartesian core: half-width and number of cells per axis."""
    h = min(p["h_iface"], p["h_wake_cross"], p["h_wake_near"])
    cover = p["r_iface_out"] + p["iface_margin"]
    if p["fine_half"]:
        if p["fine_half"] < cover - 1e-12:
            raise ValueError(f"fine_half = {p['fine_half']:g} does not cover the interface shell: it must be "
                             f">= r_iface_out + iface_margin = {cover:g}")
        half = p["fine_half"]
        n = math.ceil(2 * half / h - 1e-9)
    else:
        n = 2 * math.ceil(cover / h - 1e-9)
        half = n * h / 2
    if (n + 2) ** 3 > p["max_elements"]:  # the fine cube and one side element each way
        raise ValueError(f"the layout has more than {(n + 2) ** 3} elements (max_elements = {p['max_elements']}): "
                         "check the parameters, or raise max_elements if this is intended")
    return {"half": half, "n": n, "h": 2 * half / n}


# O-grid levels. A level is a closed surface with n x n cells on each of its
# 6 faces, corners (+-h, +-h, +-h) and face-centre radius f (h <= f <=
# sqrt(3) h). f = h is a flat cube; otherwise each face is a spherical cap
# through its 4 corners, centred on the axis at distance D behind the origin,
# and each edge is the circle where two neighbouring caps meet (it lies in a
# diagonal plane through the origin, so the separators between the 6 blocks
# of a layer are planar). Corners are indexed by sign triples s (corner
# (s_x h, s_y h, s_z h)), edges join triples that differ in one sign, and
# faces are indexed by (axis, sign).
SIGNS = [(sx, sy, sz) for sx in (-1, 1) for sy in (-1, 1) for sz in (-1, 1)]
CUBE_EDGES = [(s, t) for s in SIGNS for t in SIGNS if s < t and sum(u != v for u, v in zip(s, t)) == 1]
FACE_ORDER = [(i, sign) for i in range(3) for sign in (-1, 1)]
# Nodes of a hexahedron of a layer as (face row, face column, level) offsets
# on the half-step grids of level_nodes(..., 2 n): corners and mid-edge nodes
# in gmsh order, then face and body centres.
_CORNERS = 2 * np.array([(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)])
_HEX20 = np.vstack([_CORNERS, (_CORNERS[check_mesh.EDGES[:, 0]] + _CORNERS[check_mesh.EDGES[:, 1]]) // 2])
HEX27 = np.vstack([_HEX20, [o for o in itertools.product(range(3), repeat=3) if list(o) not in _HEX20.tolist()]])


def cap_offset(f, h):
    """D of the level (f, h): its +x face is the sphere centred at (-D, 0, 0) through the corners (0 if flat)."""
    return 0.0 if f <= h * (1 + 1e-9) else (3 * h ** 2 - f ** 2) / (2 * (f - h))


def level_nodes(f, h, m):
    """Nodes of the level (f, h) on an (m + 1) x (m + 1) grid on each face, as
    gmsh meshes the transfinite surfaces of add_level (verified by
    check_nodes): the Coons patch of the face's 4 edge arcs (nodes uniform in
    angle) at parameters (i/m, j/m), projected radially onto the face's
    sphere. Returns an array (6, m + 1, m + 1, 3), faces in FACE_ORDER."""
    D = cap_offset(f, h)
    t = np.linspace(0, 1, m + 1)
    corner = lambda sy, sz: np.array([h, sy * h, sz * h])  # noqa: E731  (corners of the +x face)

    def edge(s0, s1, centre):
        p, q = corner(*s0), corner(*s1)
        if not D:
            return p + np.outer(t, q - p)
        a, b = p - np.array(centre), q - np.array(centre)
        w = b - a * (a @ b) / (a @ a)
        w *= np.linalg.norm(a) / np.linalg.norm(w)
        th = math.atan2(np.linalg.norm(np.cross(a, b)), a @ b)
        return centre + np.outer(np.cos(th * t), a) + np.outer(np.sin(th * t), w)

    z0 = edge((-1, -1), (1, -1), [-D / 2, 0, D / 2])  # edges at z = -h, +h (along y) and y = -h, +h
    z1 = edge((-1, 1), (1, 1), [-D / 2, 0, -D / 2])
    y0 = edge((-1, -1), (-1, 1), [-D / 2, D / 2, 0])
    y1 = edge((1, -1), (1, 1), [-D / 2, -D / 2, 0])
    u, v = t[:, None, None], t[None, :, None]
    S = ((1 - v) * z0[:, None] + v * z1[:, None] + (1 - u) * y0[None] + u * y1[None]
         - ((1 - u) * (1 - v) * corner(-1, -1) + u * (1 - v) * corner(1, -1)
            + (1 - u) * v * corner(-1, 1) + u * v * corner(1, 1)))
    if D:
        c = np.array([-D, 0.0, 0.0])
        S = c + (S - c) * ((f + D) / np.linalg.norm(S - c, axis=2))[..., None]
    faces = np.zeros((6,) + S.shape)
    for q, (i, sign) in enumerate(FACE_ORDER):
        faces[q][..., i], faces[q][..., (i + 1) % 3], faces[q][..., (i + 2) % 3] = sign * S[..., 0], S[..., 1], S[..., 2]
    return faces


def layer_elements(L0, L1):
    """Nodes (6 n^2, 27, 3) of the one-cell-thick layer of hexahedra between
    two levels given by level_nodes(..., 2 n); the radial edges are straight."""
    n = (L0.shape[1] - 1) // 2
    S = np.stack([L0, 0.5 * (L0 + L1), L1])
    i = 2 * np.arange(n)
    P = S[HEX27[:, 2], np.arange(6)[:, None, None, None], i[:, None, None] + HEX27[:, 0],
          i[None, :, None] + HEX27[:, 1]]
    return P.reshape(-1, 27, 3)


def zone_excess(P, t):
    """Largest ratio of an edge (or edge aspect) to its target over the
    interface-zone and wake-zone elements among the hexahedra P (check_mesh
    definitions); <= 1 meets the targets."""
    elen, cross = check_mesh.element_edges(P)
    shell, wake = check_mesh.zones(P, t)
    ex = [0.0]
    if shell.any():
        ex += [elen[shell].max() / t["h_iface"], (elen.max(axis=1) / elen.min(axis=1))[shell].max() / t["aspect_iface"]]
    if wake.any():
        stream = np.where(cross, 0.0, elen).max(axis=1) / check_mesh.wake_envelope(P[:, :, 1].max(axis=1), t)
        ex += [stream[wake].max(), np.where(cross, elen, 0.0).max(axis=1)[wake].max() / t["h_wake_cross"]]
    return max(ex)


def largest(ok, lo, hi, steps=40):
    """Largest x in [lo, hi] with ok(x) (bisection; ok is true below some
    threshold and false above it); None if ok(lo) is false."""
    if ok(hi):
        return hi
    if not ok(lo):
        return None
    for _ in range(steps):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if ok(mid) else (lo, mid)
    return lo


def ogrid_candidate(p, n, frac):
    """The O-grid core with n cells per patch edge and core-box half-width
    a = frac n h_wake_cross / 2, or a string saying why there is none.

    Level 1 has face-centre radius f1 = r_iface_out + iface_margin; the inner
    cube half-width c = f1 n / (n + 2) makes the cube cells and the first
    layer at the face centres both t1 = 2c/n thick. Of the first levels whose
    layer meets the interface targets the flattest is used (largest corner
    half-width h1 <= f1), and its corner thickness h1 - c must be at least
    t1 / growth (a strongly bulged level 1 gives thin, sheared corner cells).

    Outer levels: the corner half-widths follow h_k = h1 + lam (f_k - f1),
    lam = (a - h1) / (a - f1), so the bulge f_k - h_k shrinks linearly to 0
    at the box and the corner thickness is lam times the face-centre one.
    The face-centre thicknesses grow from t1 by <= growth (graded_sizes; the
    first one also by <= growth from the layer-1 corner thickness) under the
    largest thickness whose wake-zone edges meet the wake targets, tabulated
    at WAKE_TABLE + 1 radii from the exact level geometry. A level that admits
    no such layer (its cells are already as long as the wake limit, which
    happens near the box when frac is close to 1) rules the candidate out.

    The candidate is also ruled out if its box does not fit the domain, if an
    edge of the O-grid (exact level geometry) or of the tensor-grid axes
    around it is below edge_min_fraction h_iface, or if neighbouring sizes
    along an axis differ by more than growth (e.g. h_core against the first
    side element, which starts from the last layer's thickness).
    """
    g, f1, t = p["growth"], p["r_iface_out"] + p["iface_margin"], targets(p)
    e_min = p["edge_min_fraction"] * p["h_iface"] - check_mesh.TOL
    ok_layer = lambda L0, L1: zone_excess(layer_elements(L0, L1), t) <= SIZE_MARGIN  # noqa: E731
    c = f1 * n / (n + 2)
    t1, a = 2 * c / n, frac * n * p["h_wake_cross"] / 2
    if t1 > p["h_iface"]:
        return "cube cells coarser than h_iface"
    if a < f1 + t1:
        return "core box too small"
    if a >= min(p["W"] / 2, p["Lu"], p["Ld"]):
        return f"core box half-width {a:g} does not fit the domain"
    L0 = level_nodes(c, c, 2 * n)
    h1 = largest(lambda h: ok_layer(L0, level_nodes(f1, h, 2 * n)), c + t1 / g, f1)
    if h1 is None:
        return "no first level meets the interface targets"
    lam = (a - h1) / (a - f1)
    s_tab = np.linspace(0, a - f1, WAKE_TABLE + 1)
    t_max = []
    for s in s_tab:
        f, h = f1 + s, h1 + lam * s
        Ls = level_nodes(f, h, 2 * n)
        d = largest(lambda d: ok_layer(Ls, level_nodes(f + d, h + lam * d, 2 * n)), 0.05 * t1, 2 * a / n, 30)
        if d is None:
            return f"no layer meets the wake targets at radius {f:.3f}"
        t_max.append(d)
    t_max = np.minimum.accumulate(np.array(t_max)[::-1])[::-1]  # non-decreasing lower bound
    ds = (a - f1) / WAKE_TABLE
    try:
        outer = graded_sizes(a - f1, min(t1, (h1 - c) / lam), lambda s: t_max[min(int(s / ds), WAKE_TABLE)], g,
                             MAX_OGRID_LAYERS)
    except ValueError:
        return f"more than {MAX_OGRID_LAYERS} outer layers"
    f = np.concatenate([[c, f1], f1 + np.cumsum(outer)])
    h = np.concatenate([[c, h1], h1 + lam * (f[2:] - f1)])
    f[-1] = h[-1] = a
    levels = [level_nodes(fk, hk, 2 * n) for fk, hk in zip(f, h)]
    edge = min(float(check_mesh.element_edges(layer_elements(L0, L1))[0].min()) for L0, L1 in zip(levels, levels[1:]))
    if edge < e_min:
        return f"O-grid edges of {edge:.4f} < edge_min_fraction h_iface"
    core = {"n": n, "cube_half": c, "h_cube": t1, "a": a, "box_fraction": frac, "h_core": 2 * a / n,
            "corner_slope": lam, "levels": [(float(fk), float(hk)) for fk, hk in zip(f, h)],
            "layer_thickness_face": np.diff(f), "layer_thickness_corner": np.diff(h), "edge_min": edge,
            "wake_limit": {"s": s_tab, "thickness": t_max}}
    try:
        ax = axes(p, core)
    except ValueError as e:
        return str(e)
    side = min(float(axis.sizes.min()) for axis in ax)
    if side < e_min:
        return f"tensor-grid elements of {side:.4f} < edge_min_fraction h_iface (the domain is too small)"
    ratio = max(axis.ratio() for axis in ax)
    if ratio > g + check_mesh.TOL:
        return f"neighbouring sizes along an axis differ by {ratio:.3f}x > growth"
    core["n_elements"] = n_elements(ax, core, p)
    return core


def ogrid_layout(p):
    """The O-grid core (see ogrid_candidate) with the fewest elements in the
    whole mesh, over the core-box fractions BOX_FRACTIONS and the first
    N_OGRID_TRIES cell counts n whose cube cells meet h_iface and whose
    largest core box holds level 1 and one layer; no candidate is evaluated
    once the element count must exceed max_elements."""
    f1 = p["r_iface_out"] + p["iface_margin"]
    n = max(2, math.ceil(2 * f1 / p["h_iface"] - 2 - 1e-9))  # smallest n with t1 = 2 f1 / (n + 2) <= h_iface
    best, rejected, tried = None, [], 0
    while tried < N_OGRID_TRIES:
        if (n + 2) ** 3 + 12 * n ** 2 > p["max_elements"]:  # one side element each way, and two layers
            rejected.append(f"n >= {n}: more than max_elements = {p['max_elements']} elements (raise max_elements "
                            "if this is intended)")
            break
        if max(BOX_FRACTIONS) * n * p["h_wake_cross"] / 2 >= f1 + 2 * f1 / (n + 2):  # a >= f1 + t1
            tried += 1
            for frac in BOX_FRACTIONS:
                core = ogrid_candidate(p, n, frac)
                if isinstance(core, str):
                    rejected.append(f"n = {n}, box fraction {frac}: {core}")
                elif best is None or core["n_elements"] < best["n_elements"]:
                    best = core
        n += 1
    if best is None:
        raise ValueError("no O-grid layout meets the targets (change h_iface, h_wake_cross, aspect_iface or "
                         "edge_min_fraction, enlarge the domain, or use core = cartesian):\n  " + "\n  ".join(rejected))
    best["rejected"] = rejected
    return best


def ogrid_nodes(core):
    """All Hex8 mesh nodes of the O-grid inside the core box: the inner cube's
    uniform grid and the nodes of levels 1 .. K-1."""
    n, c = core["n"], core["cube_half"]
    u = np.linspace(-c, c, n + 1)
    cube = np.stack(np.meshgrid(u, u, u, indexing="ij"), axis=3).reshape(-1, 3)
    return np.vstack([cube] + [level_nodes(f, h, n).reshape(-1, 3) for f, h in core["levels"][1:-1]])


def axes(p, core):
    """The x, y and z axes of the tensor-product mesh (z is the same as x);
    ValueError if they need more than max_elements elements. Beside the core
    block the x (and z) sizes stay <= h_wake_cross out to r_wake."""
    if p["core"] == "cartesian":
        half, n, h_prev = core["half"], core["n"], core["h"]
    else:
        half, n, h_prev = core["a"], core["n"], core["layer_thickness_face"][-1]
    far = lambda s: p["h_far"]  # noqa: E731
    side = lambda s: p["h_wake_cross"] if half + s <= p["r_wake"] + 1e-6 else p["h_far"]  # noqa: E731
    x = Axis(-p["W"] / 2, p["W"] / 2, half, n, h_prev, side, side, p["growth"],
             math.isqrt(p["max_elements"] // (n + 2)))  # y has >= n + 2 elements
    y = Axis(-p["Ld"], p["Lu"], half, n, h_prev, wake_target(p, -half), far, p["growth"],
             p["max_elements"] // len(x.sizes) ** 2)
    return x, y, x


def block_grid(ax, hole=None):
    """Build the transfinite block grid spanned by the breakpoints of the three
    axes (gmsh built-in kernel), leaving out the block `hole` (index triple).
    Returns the volume tags, for each axis d and side (0 = low, 1 = high) the
    tags of the boundary surfaces, and the point, line and surface tags."""
    g = gmsh.model.geo
    nb = tuple(len(a.breaks) for a in ax)
    unit = np.eye(3, dtype=int)
    P = np.zeros(nb, dtype=int)
    for idx in np.ndindex(nb):
        P[idx] = g.addPoint(*(ax[d].breaks[idx[d]] for d in range(3)))

    # Lines along axis d, from P[idx] to P[idx + e_d] (increasing coordinate).
    L = [{}, {}, {}]
    for d in range(3):
        for idx in np.ndindex(nb):
            if idx[d] < nb[d] - 1:
                tag = g.addLine(P[idx], P[tuple(idx + unit[d])])
                n_nodes, coef = ax[d].segment(idx[d])
                g.mesh.setTransfiniteCurve(tag, n_nodes, "Progression", coef)
                L[d][idx] = tag

    # Faces normal to axis d, spanned by the two other axes a and b.
    S = [{}, {}, {}]
    for d in range(3):
        a, b = (d + 1) % 3, (d + 2) % 3
        for idx in np.ndindex(nb):
            if idx[a] < nb[a] - 1 and idx[b] < nb[b] - 1:
                ia, ib = tuple(idx + unit[a]), tuple(idx + unit[b])
                loop = g.addCurveLoop([L[a][idx], L[b][ia], -L[a][ib], -L[b][idx]])
                tag = g.addPlaneSurface([loop])
                g.mesh.setTransfiniteSurface(tag)
                g.mesh.setRecombine(2, tag)
                S[d][idx] = tag

    volumes = []
    for idx in np.ndindex(tuple(n - 1 for n in nb)):
        if idx == hole:
            continue
        faces = [S[d][i] for d in range(3) for i in (idx, tuple(idx + unit[d]))]
        tag = g.addVolume([g.addSurfaceLoop(faces)])
        g.mesh.setTransfiniteVolume(tag)
        volumes.append(tag)

    boundary = [[[t for idx, t in S[d].items() if idx[d] == side * (nb[d] - 1)] for side in (0, 1)]
                for d in range(3)]
    return volumes, boundary, (P, L, S)


def face_corners(i, sign):
    """The 4 corner sign triples of the level face normal to axis i, in cyclic order."""
    j, k = (i + 1) % 3, (i + 2) % 3
    out = []
    for u, v in ((-1, -1), (1, -1), (1, 1), (-1, 1)):
        s = [0, 0, 0]
        s[i], s[j], s[k] = sign, u, v
        out.append(tuple(s))
    return out


class Curves:
    """Directed curves between point tags, for building curve loops."""

    def __init__(self):
        self.tags = {}

    def add(self, p, q, tag):
        """Record curve `tag` as running from point p to point q."""
        self.tags[p, q] = tag
        return tag

    def signed(self, p, q):
        """Signed tag of the curve p -> q (negative if it was created as q -> p)."""
        return self.tags[p, q] if (p, q) in self.tags else -self.tags[q, p]

    def quad(self, pts, sphere_center=None):
        """Transfinite, recombined 4-sided surface through the point cycle: a
        plane, or a spherical patch (surface filling on the sphere through
        its corners centred at point `sphere_center`)."""
        g = gmsh.model.geo
        loop = g.addCurveLoop([self.signed(pts[m], pts[(m + 1) % 4]) for m in range(4)])
        s = g.addPlaneSurface([loop]) if sphere_center is None else g.addSurfaceFilling([loop], sphereCenterTag=sphere_center)
        g.mesh.setTransfiniteSurface(s)
        g.mesh.setRecombine(2, s)
        return s


def transfinite_volume(faces):
    """Transfinite hex block bounded by 6 transfinite surfaces."""
    g = gmsh.model.geo
    v = g.addVolume([g.addSurfaceLoop(faces)])
    g.mesh.setTransfiniteVolume(v)
    return v


def add_level(cv, f, h, n, origin):
    """One closed O-grid level (f, h) with n x n cells per face (see
    level_nodes). Returns the corner point tags by sign triple and the face
    tags by (axis, sign)."""
    g = gmsh.model.geo
    pts = {s: g.addPoint(h * s[0], h * s[1], h * s[2]) for s in SIGNS}
    D = cap_offset(f, h)
    centres = {}

    def centre(vec):
        key = tuple(round(v, 12) for v in vec)
        if key not in centres:
            centres[key] = origin if not any(key) else g.addPoint(*vec)
        return centres[key]

    for s, t in CUBE_EDGES:
        if not D:
            tag = g.addLine(pts[s], pts[t])
        else:  # the edge is shared by the caps normal to the two axes where s == t
            tag = g.addCircleArc(pts[s], centre([-0.5 * D * s[i] if s[i] == t[i] else 0.0 for i in range(3)]), pts[t])
        g.mesh.setTransfiniteCurve(tag, n + 1)
        cv.add(pts[s], pts[t], tag)
    faces = {}
    for i, sign in FACE_ORDER:
        ctr = centre([-D * sign if j == i else 0.0 for j in range(3)]) if D else None
        faces[i, sign] = cv.quad([pts[s] for s in face_corners(i, sign)], ctr)
    return pts, faces


def box_level(cv, grid, hole):
    """The faces of block `hole` of the tensor grid as the outermost O-grid level."""
    P, L, S = grid
    corner = lambda s: tuple(hole[d] + (s[d] > 0) for d in range(3))  # noqa: E731
    pts = {s: int(P[corner(s)]) for s in SIGNS}
    for s, t in CUBE_EDGES:
        d = next(i for i in range(3) if s[i] != t[i])
        lo, hi = (s, t) if s[d] < t[d] else (t, s)
        cv.add(pts[lo], pts[hi], L[d][corner(lo)])
    faces = {}
    for i, sign in FACE_ORDER:
        idx = list(hole)
        idx[i] += sign > 0
        faces[i, sign] = S[i][tuple(idx)]
    return pts, faces


def build_ogrid(core, grid, hole):
    """Fill block `hole` of the tensor grid with the O-grid: the inner cube and
    one layer of 6 transfinite blocks between consecutive levels."""
    g = gmsh.model.geo
    cv = Curves()
    origin = g.addPoint(0, 0, 0)
    n = core["n"]
    levels = [add_level(cv, f, h, n, origin) for f, h in core["levels"][:-1]] + [box_level(cv, grid, hole)]
    volumes = [transfinite_volume(list(levels[0][1].values()))]
    for (p0, f0), (p1, f1) in zip(levels, levels[1:]):
        for s in SIGNS:  # corner rays, one element each
            g.mesh.setTransfiniteCurve(cv.add(p0[s], p1[s], g.addLine(p0[s], p1[s])), 2)
        sep = {(s, t): cv.quad([p0[s], p0[t], p1[t], p1[s]]) for s, t in CUBE_EDGES}
        for i, sign in FACE_ORDER:
            cc = face_corners(i, sign)
            side = [sep[tuple(sorted((cc[m], cc[(m + 1) % 4])))] for m in range(4)]
            volumes.append(transfinite_volume([f0[i, sign], f1[i, sign]] + side))
    return volumes


def n_elements(ax, core, p):
    """Designed element count (the O-grid replaces the core block's n^3)."""
    n = int(np.prod([len(a.sizes) for a in ax]))
    if p["core"] == "ogrid":
        n += 6 * core["n"] ** 2 * (len(core["levels"]) - 1)
    return n


def nearest_distance(A, B, chunk=512):
    """Distance from each point of A to the nearest point of B."""
    b2 = (B ** 2).sum(axis=1)
    out = [np.sqrt(np.maximum(((a ** 2).sum(axis=1)[:, None] + b2 - 2 * a @ B.T).min(axis=1), 0))
           for a in np.array_split(A, max(1, len(A) // chunk))]
    return np.concatenate(out)


def check_nodes(ax, core, p):
    """Verify that gmsh produced only hexahedra, as many as designed, and that
    every mesh node lies on the design: on the tensor-product distributions
    outside the O-grid, on the level nodes of the model inside it."""
    types, tags, _ = gmsh.model.mesh.getElements(3)
    if list(types) != [5]:
        raise RuntimeError(f"expected only 8-node hexahedra, got element types {list(types)}")
    n_expected = n_elements(ax, core, p)
    if len(tags[0]) != n_expected:
        raise RuntimeError(f"{len(tags[0])} hexahedra, expected {n_expected}")
    _, xyz, _ = gmsh.model.mesh.getNodes(3, -1, includeBoundary=True)
    xyz = xyz.reshape(-1, 3)
    if p["core"] == "ogrid":
        inner = np.abs(xyz).max(axis=1) < core["a"] * (1 - 1e-9)
        err = nearest_distance(xyz[inner], ogrid_nodes(core)).max()
        if err > 1e-6:
            raise RuntimeError(f"O-grid mesh nodes deviate from the design by {err:.2e}")
        xyz = xyz[~inner]
    for d in range(3):
        err = np.abs(xyz[:, d, None] - ax[d].nodes[None, :]).min(axis=1).max()
        if err > 1e-6:
            raise RuntimeError(f"axis {d}: mesh nodes deviate from the design by {err:.2e}")


def build(p):
    """Mesh the domain, write p['name'].msh and return the core layout, the
    axes and the check_mesh metrics."""
    validate(p)
    core = cartesian_layout(p) if p["core"] == "cartesian" else ogrid_layout(p)
    ax = axes(p, core)
    n_design = n_elements(ax, core, p)
    if n_design > p["max_elements"]:
        raise ValueError(f"the layout has {n_design} elements (max_elements = {p['max_elements']}); "
                         "check the parameters, or raise max_elements if this is intended")
    gmsh.initialize()
    try:
        gmsh.option.setNumber("General.Terminal", 1)
        gmsh.option.setNumber("General.Verbosity", 2)  # warnings and errors only
        gmsh.option.setNumber("General.NumThreads", p["threads"])
        gmsh.model.add(p["name"])
        hole = tuple(a.core for a in ax) if p["core"] == "ogrid" else None
        volumes, bnd, grid = block_grid(ax, hole)
        if hole:
            volumes += build_ogrid(core, grid, hole)
        gmsh.model.geo.synchronize()

        # Order matters: gmsh2nek assigns boundary IDs 1..6 in this order.
        for name, d, side in (("inflow", 1, 1), ("outflow", 1, 0), ("xmin", 0, 0),
                              ("xmax", 0, 1), ("zmin", 2, 0), ("zmax", 2, 1)):
            gmsh.model.addPhysicalGroup(2, bnd[d][side], name=name)
        gmsh.model.addPhysicalGroup(3, volumes, name="fluid")

        gmsh.model.mesh.generate(3)
        check_nodes(ax, core, p)
        gmsh.model.mesh.setOrder(2)
        metrics = check_mesh.measure(check_mesh.read_model(), targets(p))
        gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
        gmsh.option.setNumber("Mesh.Binary", 1)
        gmsh.write(p["name"] + ".msh")
    finally:
        gmsh.finalize()
    return core, ax, metrics


# Which parameters to change when a target fails; the other checks (element
# type, conformity, groups, periodicity) failing means a construction bug.
ADVICE = {
    "negative_jacobians": "growth (smaller)",
    "minSJ": "growth (smaller), or core = cartesian",
    "interface_edge_max": "h_iface",
    "interface_aspect_max": "aspect_iface or h_iface",
    "wake_cross_max": "h_wake_cross or r_wake",
    "wake_stream_envelope_excess": "h_wake_near, h_wake_far or y_wake_near / y_wake_ramp",
    "wake_near_stream_max": "h_wake_near",
    "wake_far_stream_max": "h_wake_far",
    "far_edge_max": "h_far",
    "hJ_ratio_max": "growth (smaller)",
    "thickness_ratio_max": "growth (smaller)",
    "edge_min": "edge_min_fraction, or the domain (W, Lu, Ld) and h_far against h_iface",
    "axis_ratio_max": "the domain (W, Lu, Ld): the outermost elements are squeezed",
}


def verify(res, name):
    """Raise if a check failed, naming the targets and the parameters to change."""
    failed = [k for k, (_, _, ok) in res.items() if not ok]
    if failed:
        lines = [f"  {k} = {res[k][0]} (limit {res[k][1]}): change {ADVICE.get(k, 'nothing; this is a bug')}"
                 for k in failed]
        raise RuntimeError(f"{name}.msh misses {len(failed)} target(s) (written for inspection; "
                           f"see {name}.plan.json):\n" + "\n".join(lines))


def json_default(o):
    """JSON encoding of numpy arrays and scalars."""
    return o.tolist() if hasattr(o, "tolist") else str(o)


def generate(p, preset="custom"):
    """Mesh with the parameters p, check the result and write NAME.msh and
    NAME.plan.json; raise if a target fails. Returns the plan."""
    t0 = time.time()
    core, ax, metrics = build(p)
    res = check_mesh.evaluate(metrics, targets(p))
    ratio = max(a.ratio() for a in ax)
    res["axis_ratio_max"] = (ratio, p["growth"], ratio <= p["growth"] + check_mesh.TOL)
    ne = metrics["n_elements"]
    plan = {"preset": preset, "parameters": p,
            "layout": {"core": core, "x_and_z": ax[0].plan(), "y": ax[1].plan(),
                       "blocks": [len(a.breaks) - 1 for a in ax]},
            "n_elements": ne,
            "n_elements_core": ne - int(np.prod([len(a.sizes) for a in ax])) + core["n"] ** 3,
            "gll_points": {f"N{N}": check_mesh.gll_points(ne, N) for N in sorted({5, 7, 9, p["N_target"]})},
            "nekrs_ls": metrics.get("nekrs_ls"),
            "checks": {k: {"value": v, "limit": lim, "ok": ok} for k, (v, lim, ok) in res.items()},
            "generation_time_s": round(time.time() - t0, 2),
            "metrics": metrics}
    with open(p["name"] + ".plan.json", "w") as f:
        json.dump(plan, f, indent=1, default=json_default)

    np.set_printoptions(precision=3, linewidth=150)
    if p["core"] == "ogrid":
        print(f"O-grid: n = {core['n']}, inner cube half-width {core['cube_half']:.4f} (cells {core['h_cube']:.4f}), "
              f"box half-width {core['a']:g}, {len(core['levels']) - 1} layers")
        print("  levels (face-centre radius, corner half-width):", np.array(core["levels"]).round(4).tolist())
    else:
        print(f"fine cube: half-width {core['half']:g}, {core['n']} cells of {core['h']:.4g}")
    for label, a in (("x = z", ax[0]), ("y", ax[1])):
        print(f"{label}: {len(a.sizes)} elements in {len(a.breaks) - 1} segments\n  sizes: {a.sizes}")
    print(f"{p['name']}.msh ({plan['generation_time_s']} s)\n{check_mesh.summary(metrics)}\n{check_mesh.table(res)}")
    verify(res, p["name"])
    return plan


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0], allow_abbrev=False,
                                 epilog="presets: " + "; ".join(f"{k} {v}" for k, v in PRESETS.items()))
    ap.add_argument("--preset", choices=list(PRESETS), default="default")
    for k, v in PARAMS.items():
        ap.add_argument(f"--{k}", type=type(v), choices=CORES if k == "core" else None, help=f"default {v}")
    ap.add_argument("--print_name", action="store_true", help="print NAME and exit (used by ./mesh)")
    a = ap.parse_args()
    p = {**PARAMS, **PRESETS[a.preset], **{k: getattr(a, k) for k in PARAMS if getattr(a, k) is not None}}
    if a.print_name:
        print(p["name"])
        return
    check_mesh.limit_memory(p["max_memory_gb"])
    generate(p, a.preset)


if __name__ == "__main__":
    main()
