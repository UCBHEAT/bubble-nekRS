"""Resolution, quality and boundary checks of the 3D rising-bubble hex mesh.

usage: check_mesh.py NAME.msh [--re2 NAME.re2] [--json OUT.json]
                     [--plan NAME.plan.json] [--W 6 --h_iface 0.2 ...]

Conventions: y is vertical (gravity -y), the bubble (D = 1) is centred at the
origin, the wake is below it (y < 0), and x, z are periodic with period W.
Elements are classified from all their nodes:

- interface zone: elements overlapping the shell r_iface_in <= r <= r_iface_out
  (r = |x|)
- wake zone: other elements overlapping the column rho <= r_wake
  (rho = sqrt(x^2 + z^2)), y_wake_end <= y <= y_wake_top
- far field: everything else

Edge lengths follow the curved edge (|a - mid| + |mid - b|) of second-order
elements. An edge is streamwise if its chord is within 60 degrees of the y
axis, else cross-stream. "h_J" is the element volume^(1/3) (from the
Jacobian), the length scale nekRS-LS uses for its interface width eps. The
wake targets are checked in bands (an element counts in every band its
y-extent reaches) and against the linear envelope of the streamwise target
evaluated at the top of each element. For information it also reports the
classical corner scaled Jacobian and a vertical-flow CFL length per element
(both are small in sheared cells, e.g. at O-grid corners).

The targets default to TARGETS below; --plan reads them from the plan file
written by generate_mesh.py, and --KEY VALUE overrides single ones. With
--re2 the gmsh2nek output is read and checked against the .msh as well.
Prints a summary and a table of the targets (--json writes all metrics);
exit status 1 if any target fails.
"""
import argparse
import json
import math
import resource
import sys

import gmsh
import numpy as np

# Target values (lengths in bubble diameters); generate_mesh.py uses the same
# names for its parameters.
TARGETS = {
    # Domain: x, z in [-W/2, W/2] (periodic), y in [-Ld, +Lu].
    "W": 6.0,
    "Lu": 5.0,
    "Ld": 10.0,
    # Interface zone: largest edge and largest edge aspect ratio.
    "r_iface_in": 0.3,
    "r_iface_out": 0.75,
    "h_iface": 0.2,
    "aspect_iface": 2.0,
    # Wake zone: cross-stream edges <= h_wake_cross; streamwise <= h_wake_near
    # down to y_wake_near, a linear ramp to h_wake_far at y_wake_ramp, then
    # h_wake_far down to y_wake_end.
    "r_wake": 0.75,
    "y_wake_top": -0.5,
    "h_wake_cross": 0.3,
    "h_wake_near": 0.25,
    "y_wake_near": -1.5,
    "h_wake_far": 0.5,
    "y_wake_ramp": -5.0,
    "y_wake_end": -6.0,
    # Far field: largest edge.
    "h_far": 1.25,
    # Grading and quality: face-neighbour h_J ratio, face-normal thickness
    # ratio, gmsh minimum scaled Jacobian, and the smallest edge anywhere as a
    # fraction of h_iface ("no tiny cells").
    "hJ_ratio_max": 1.5,
    "thickness_ratio_max": 2.0,
    "minSJ_min": 0.5,
    "edge_min_fraction": 0.5,
    # nekRS-LS settings are derived for N = 7 and N_target (Weber number and
    # liquid/gas density ratio of the case).
    "N_target": 7,
    "We": 2.667,
    "rho_ratio": 40.0,
}
GROUPS = ("inflow", "outflow", "xmin", "xmax", "zmin", "zmax")  # gmsh2nek IDs 1..6
TOL = 1e-9  # slack when comparing lengths with targets

# gmsh hex numbering: edges as corner pairs (mid-edge node 8 + m of edge m),
# faces as corner quadruples, and the opposite of each face.
EDGES = np.array([(0, 1), (0, 3), (0, 4), (1, 2), (1, 5), (2, 3), (2, 6), (3, 7), (4, 5), (4, 7), (5, 6), (6, 7)])
FACES = np.array([(0, 3, 2, 1), (0, 1, 5, 4), (0, 4, 7, 3), (1, 2, 6, 5), (2, 3, 7, 6), (4, 5, 6, 7)])
OPPOSITE = np.array([5, 4, 3, 2, 1, 0])
# The neighbours of each corner along its three edges, as a right-handed triple.
CORNER_NEIGHBOURS = np.array([(1, 3, 4), (2, 0, 5), (3, 1, 6), (0, 2, 7), (7, 5, 0), (4, 6, 1), (5, 7, 2), (6, 4, 3)])


def read_model():
    """Volume elements, nodes and boundary groups of the current gmsh model."""
    ntags, xyz, _ = gmsh.model.mesh.getNodes()
    X = np.zeros((int(ntags.max()) + 1, 3))
    X[ntags.astype(np.int64)] = xyz.reshape(-1, 3)
    types = gmsh.model.mesh.getElementTypes(3)
    counts = {gmsh.model.mesh.getElementProperties(t)[0]: len(gmsh.model.mesh.getElementsByType(t)[0])
              for t in types}
    hexes = [t for t in types if gmsh.model.mesh.getElementProperties(t)[0].startswith("Hexahedron")]
    if not hexes:
        raise RuntimeError(f"no hexahedra in the mesh (element types {counts})")
    et = max(hexes, key=lambda t: counts[gmsh.model.mesh.getElementProperties(t)[0]])
    etags, enodes = gmsh.model.mesh.getElementsByType(et)
    nper = gmsh.model.mesh.getElementProperties(et)[3]
    loc, w = gmsh.model.mesh.getIntegrationPoints(et, "Gauss4")
    _, det, _ = gmsh.model.mesh.getJacobians(et, loc)
    groups = {}
    for dim, tag in gmsh.model.getPhysicalGroups(2):
        faces = []
        for ent in gmsh.model.getEntitiesForPhysicalGroup(dim, tag):
            for ft, _, fn in zip(*gmsh.model.mesh.getElements(2, ent)):
                k = gmsh.model.mesh.getElementProperties(ft)[3]
                faces.append(X[fn.astype(np.int64).reshape(-1, k)[:, :4]].mean(axis=1))
        groups[gmsh.model.getPhysicalName(dim, tag)] = np.vstack(faces) if faces else np.zeros((0, 3))
    return {"X": X, "conn": enodes.reshape(len(etags), nper).astype(np.int64), "element_types": counts,
            "det": det.reshape(len(etags), -1), "weights": np.asarray(w),
            "minSJ": np.array(gmsh.model.mesh.getElementQualities(etags, "minSJ")),
            "minSICN": np.array(gmsh.model.mesh.getElementQualities(etags, "minSICN")),
            "groups": groups,
            "volume_groups": [gmsh.model.getPhysicalName(3, t) for _, t in gmsh.model.getPhysicalGroups(3)]}


def load_msh(path):
    """read_model() of a .msh file."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        gmsh.open(path)
        return read_model()
    finally:
        gmsh.finalize()


def _stats(v):
    return {"min": float(v.min()), "p01": float(np.percentile(v, 1)), "median": float(np.median(v)),
            "max": float(v.max())}


def _where(P, k):
    """Centroid and radius of element k, for reports."""
    c = P[k, :8].mean(axis=0)
    return {"centroid": [round(float(v), 3) for v in c], "r": round(float(np.linalg.norm(c)), 3)}


def gll_points(n_elements, N):
    return int(n_elements * (N + 1) ** 3)


def nekrs_ls(iz, hJ, cfl_length, N, t):
    """nekRS-LS settings for polynomial order N from the interface-zone
    metrics iz and the h_J (element volume^(1/3)) and CFL lengths of all
    elements: interface width eps = 1.5 h_J / N from the median
    interface-zone h_J, and for comparison nekRS-LS's default
    (interfaceWidthFactor 1: the largest 2 J^(1/3) over all GLL points / N,
    which is the largest h_J / N for affine elements), the explicit
    surface-tension (capillary) dt limit at the smallest interface-zone h_J,
    the advective dt at CFL 0.5 for a unit vertical velocity, and the element
    size range (the reinit pseudo-step count scales with Hmax/Hmin)."""
    eps = 1.5 * iz["hJ_median"] / N
    sigma = 1.0 / t["We"]
    capillary = 0.8 * np.sqrt((1 + 1 / t["rho_ratio"]) * (iz["hJ_min"] / N) ** 3 / (4 * np.pi * sigma))
    return {"N": N, "gll_points": gll_points(len(hJ), N),
            "recommended_interfaceWidthValue": eps, "R_over_eps": 0.5 / eps,
            "default_eps_hJmax_over_N": hJ.max() / N, "R_over_default_eps": 0.5 / (hJ.max() / N),
            "capillary_dt_0.8x": float(capillary),
            "advective_dt_cfl0.5_U1": 0.5 * float(cfl_length.min()) * gll_spacing(N),
            "Hmax_over_Hmin": float(hJ.max() / hJ.min())}


def element_edges(P):
    """Edge lengths (ne, 12) of hexahedra given by their nodes P (ne, >= 8, 3)
    in gmsh order, along the curved edge when the mid-edge nodes (8..19) are
    given, and whether each edge is cross-stream (chord more than 60 degrees
    from the y axis)."""
    C = P[:, :8]
    chord = C[:, EDGES[:, 1]] - C[:, EDGES[:, 0]]
    if P.shape[1] >= 20:
        M = P[:, 8:20]
        elen = np.linalg.norm(M - C[:, EDGES[:, 0]], axis=2) + np.linalg.norm(C[:, EDGES[:, 1]] - M, axis=2)
    else:
        elen = np.linalg.norm(chord, axis=2)
    return elen, np.abs(chord[:, :, 1]) < 0.5 * np.linalg.norm(chord, axis=2)


def zones(P, t):
    """Masks of the interface-zone and wake-zone elements (by overlap, from
    all given nodes P (ne, k, 3)); the rest is far field."""
    r = np.linalg.norm(P, axis=2)
    rho = np.hypot(P[:, :, 0], P[:, :, 2])
    shell = (r.min(axis=1) <= t["r_iface_out"]) & (r.max(axis=1) >= t["r_iface_in"])
    wake = (~shell & (rho.min(axis=1) <= t["r_wake"]) & (P[:, :, 1].min(axis=1) <= t["y_wake_top"])
            & (P[:, :, 1].max(axis=1) >= t["y_wake_end"]))
    return shell, wake


def corner_quality(C):
    """Per element: the classical (Verdict) scaled Jacobian, i.e. the smallest
    det(e1, e2, e3) / (|e1| |e2| |e3|) over the 8 corners (gmsh's minSJ is
    normalised differently), and the vertical-flow CFL length 1 / max over
    the corners of sum_i |grad xi_i . e_y| (xi_i in [0, 1]; the edge length
    for a cube), which is short in sheared cells."""
    E = [C[:, nb] - C[:, [k]] for k, nb in enumerate(CORNER_NEIGHBOURS)]   # 8 x (ne, 3, 3)
    E = np.stack(E, axis=1)                                                 # (ne, 8, 3 edges, 3)
    e1, e2, e3 = E[:, :, 0], E[:, :, 1], E[:, :, 2]
    det = np.einsum("ijk,ijk->ij", e1, np.cross(e2, e3))
    sign = np.sign(det.sum(axis=1, keepdims=True))
    sj = (sign * det / np.prod(np.linalg.norm(E, axis=3), axis=2)).min(axis=1)
    rate = (np.abs(np.cross(e2, e3)[..., 1]) + np.abs(np.cross(e3, e1)[..., 1])
            + np.abs(np.cross(e1, e2)[..., 1])) / np.abs(det)
    return sj, 1.0 / rate.max(axis=1)


def gll_spacing(N):
    """Smallest GLL point spacing for polynomial order N, as a fraction of the element."""
    x = np.sort(np.polynomial.legendre.Legendre.basis(N).deriv().roots())
    return float((1 - x[-1]) / 2)


def measure(mesh, t):
    """Metrics of read_model() output for the targets t (JSON-ready dict)."""
    X, conn = mesh["X"], mesh["conn"]
    ne = len(conn)
    P = X[conn]                                               # (ne, nper, 3)
    C = P[:, :8]
    out = {"element_types": mesh["element_types"], "n_elements": ne}
    for N in sorted({5, 7, 9, t["N_target"]}):
        out[f"gll_points_N{N}"] = gll_points(ne, N)

    # Edges, h_J (element volume^(1/3)) and quality.
    elen, cross = element_edges(P)
    vol = mesh["det"] @ mesh["weights"]
    hJ = np.cbrt(np.abs(vol))
    aspect = elen.max(axis=1) / elen.min(axis=1)
    out["n_negative_jacobian_elements"] = int(np.sum(mesh["det"].min(axis=1) <= 0))
    out["total_volume"] = float(vol.sum())
    sj = mesh["minSJ"]
    csj, lcfl = corner_quality(C)
    out["quality"] = {"minSJ": _stats(sj), "minSICN": _stats(mesh["minSICN"]), "edge_aspect": _stats(aspect),
                      "worst_minSJ": [{"minSJ": float(sj[k]), **_where(P, k)} for k in np.argsort(sj)[:5]],
                      "corner_SJ": {**_stats(csj), "worst_at": _where(P, int(np.argmin(csj)))},
                      "cfl_length": {**_stats(lcfl), "worst_at": _where(P, int(np.argmin(lcfl)))}}
    k = int(np.argmin(elen.min(axis=1)))
    out["edge_min"] = {"value": float(elen.min()), **_where(P, k)}
    out["hJ"] = {"min": float(hJ.min()), "max": float(hJ.max()), "Hmax_over_Hmin": float(hJ.max() / hJ.min())}

    # Face neighbours: h_J ratio, and the thickness ratio normal to the shared
    # face (h_J hides one-directional jumps). An element's thickness across
    # face f is the distance between the centroids of f and its opposite face.
    fkey = np.sort(conn[:, :8][:, FACES], axis=2).reshape(-1, 4)
    _, inv, cnt = np.unique(fkey, axis=0, return_inverse=True, return_counts=True)
    order = np.argsort(inv.ravel(), kind="stable")
    first = np.nonzero(np.diff(inv.ravel()[order]) == 0)[0]
    f1, f2 = order[first], order[first + 1]                   # indices into (ne * 6)
    e1, e2 = f1 // 6, f2 // 6
    ratio = np.maximum(hJ[e1] / hJ[e2], hJ[e2] / hJ[e1])
    fc = C[:, FACES].mean(axis=2)
    thick = np.linalg.norm(fc - fc[:, OPPOSITE], axis=2).ravel()
    tr = np.maximum(thick[f1] / thick[f2], thick[f2] / thick[f1])
    k = int(np.argmax(tr))
    out["neighbor_hJ_ratio"] = {"max": float(ratio.max()), "p99": float(np.percentile(ratio, 99)),
                                "n_pairs": int(len(ratio))}
    out["neighbor_thickness_ratio"] = {
        "max": float(tr.max()), "p99": float(np.percentile(tr, 99)), "n_over_1.5": int(np.sum(tr > 1.5 * (1 + 1e-6))),
        "worst_at": [round(float(v), 3) for v in 0.5 * (C[e1[k]].mean(0) + C[e2[k]].mean(0))],
        "worst_thicknesses": [float(thick[f1[k]]), float(thick[f2[k]])]}
    out["n_faces_shared_by_more_than_2"] = int(np.sum(cnt > 2))
    out["n_boundary_faces"] = int(np.sum(cnt == 1))

    # Regions.
    shell, wake = zones(P, t)
    far = ~(shell | wake)
    ylo, yhi = P[:, :, 1].min(axis=1), P[:, :, 1].max(axis=1)
    Lc = np.where(cross, elen, 0.0).max(axis=1)               # longest cross-stream edge
    Ls = np.where(cross, 0.0, elen).max(axis=1)               # longest streamwise edge

    def worst(m, v):
        """Largest value of v over the elements m, and where it is."""
        k = int(np.argmax(np.where(m, v, -np.inf)))
        return {"value": float(v[k]), **_where(P, k)}

    def region(m):
        if not m.any():
            return None
        return {"n_elem": int(m.sum()), "edge_min": float(elen[m].min()), "edge_median": float(np.median(elen[m])),
                "edge_max": worst(m, elen.max(axis=1)), "aspect_max": worst(m, aspect),
                "cross_edge_max": worst(m, Lc), "stream_edge_max": worst(m, Ls),
                "hJ_min": float(hJ[m].min()), "hJ_median": float(np.median(hJ[m])), "hJ_max": float(hJ[m].max())}

    reg = {"interface_zone": region(shell), "wake_zone": region(wake), "far_field": region(far)}
    bands = {}
    for name, lo, hi in (("near", t["y_wake_near"], t["y_wake_top"]), ("mid1", -3.0, -2.0), ("mid2", -4.5, -3.5),
                         ("far", t["y_wake_end"], t["y_wake_ramp"])):
        m = wake & (ylo < hi) & (yhi > lo)
        if m.any():
            bands[name] = {"y": [lo, hi], "n_elem": int(m.sum()), "cross_edge_max": worst(m, Lc),
                           "stream_edge_max": worst(m, Ls)}
    reg["wake_bands"] = bands
    if wake.any():
        reg["wake_envelope_excess"] = worst(wake, Ls - wake_envelope(yhi, t))
    out["regions"] = reg

    # nekRS-LS settings.
    if reg["interface_zone"]:
        out["nekrs_ls"] = {f"N{N}": nekrs_ls(reg["interface_zone"], hJ, lcfl, N, t)
                           for N in sorted({7, t["N_target"]})}

    # Boundary groups: order, coverage, planes and periodic matching.
    g = mesh["groups"]
    out["boundary_groups"] = {name: int(len(fcs)) for name, fcs in g.items()}
    out["volume_groups"] = mesh["volume_groups"]
    plane = {"inflow": (1, t["Lu"]), "outflow": (1, -t["Ld"]), "xmin": (0, -t["W"] / 2), "xmax": (0, t["W"] / 2),
             "zmin": (2, -t["W"] / 2), "zmax": (2, t["W"] / 2)}
    out["boundary_plane_deviation"] = max((float(np.abs(g[n][:, d] - v).max()) for n, (d, v) in plane.items()
                                           if n in g and len(g[n])), default=None)
    out["periodic_match"] = {
        "x": match(g.get("xmin", np.zeros((0, 3))), g.get("xmax", np.zeros((0, 3))), (t["W"], 0, 0)),
        "z": match(g.get("zmin", np.zeros((0, 3))), g.get("zmax", np.zeros((0, 3))), (0, 0, t["W"]))}
    return out


def wake_envelope(y, t):
    """Streamwise size target at height y (array) in the wake."""
    ramp = (t["y_wake_near"] - y) / (t["y_wake_near"] - t["y_wake_ramp"])
    return t["h_wake_near"] + np.clip(ramp, 0, 1) * (t["h_wake_far"] - t["h_wake_near"])


def match(A, B, shift):
    """Largest distance between the face centres of A shifted by `shift` and
    those of B, matched by sorting on rounded coordinates."""
    if len(A) != len(B) or len(A) == 0:
        return {"n1": len(A), "n2": len(B), "max_mismatch": None}
    A2 = A + np.array(shift)
    ka = np.lexsort(np.round(A2, 6).T[::-1])
    kb = np.lexsort(np.round(B, 6).T[::-1])
    return {"n1": len(A), "n2": len(B), "max_mismatch": float(np.abs(A2[ka] - B[kb]).max())}


def evaluate(m, t):
    """Target name -> (value, limit, ok) for the metrics m of measure()."""
    res = {}

    def at_most(name, value, limit):
        if value is not None:
            res[name] = (value, limit, bool(value <= limit + TOL))

    def at_least(name, value, limit):
        res[name] = (value, limit, bool(value >= limit - TOL))

    res["hex27_only"] = (sorted(m["element_types"]), ["Hexahedron 27"], list(m["element_types"]) == ["Hexahedron 27"])
    at_most("negative_jacobians", m["n_negative_jacobian_elements"], 0)
    at_least("minSJ", m["quality"]["minSJ"]["min"], t["minSJ_min"])
    reg = m["regions"]
    iz, wz, ff, bands = reg["interface_zone"], reg["wake_zone"], reg["far_field"], reg["wake_bands"]
    if iz:
        at_most("interface_edge_max", iz["edge_max"]["value"], t["h_iface"])
        at_most("interface_aspect_max", iz["aspect_max"]["value"], t["aspect_iface"])
    if wz:
        at_most("wake_cross_max", wz["cross_edge_max"]["value"], t["h_wake_cross"])
        at_most("wake_stream_envelope_excess", reg["wake_envelope_excess"]["value"], 0.0)
    for band, limit in (("near", t["h_wake_near"]), ("far", t["h_wake_far"])):
        if band in bands:
            at_most(f"wake_{band}_stream_max", bands[band]["stream_edge_max"]["value"], limit)
    if ff:
        at_most("far_edge_max", ff["edge_max"]["value"], t["h_far"])
    at_most("hJ_ratio_max", m["neighbor_hJ_ratio"]["max"], t["hJ_ratio_max"])
    at_most("thickness_ratio_max", m["neighbor_thickness_ratio"]["max"], t["thickness_ratio_max"])
    at_least("edge_min", m["edge_min"]["value"], t["edge_min_fraction"] * t["h_iface"])
    at_most("faces_shared_by_more_than_2", m["n_faces_shared_by_more_than_2"], 0)
    names = list(m["boundary_groups"])
    res["group_order"] = (names, list(GROUPS), names == list(GROUPS) and m["volume_groups"] == ["fluid"])
    n_grouped = sum(m["boundary_groups"].values())
    res["boundary_faces_grouped"] = (n_grouped, m["n_boundary_faces"], n_grouped == m["n_boundary_faces"])
    at_most("boundary_plane_deviation", m["boundary_plane_deviation"], 1e-8)
    for d, pm in m["periodic_match"].items():
        ok = pm["max_mismatch"] is not None and pm["max_mismatch"] <= 1e-8
        res[f"periodic_{d}"] = (pm["max_mismatch"], 1e-8, ok)
    return res


def read_re2(path):
    """Header, corner Jacobians, curved edges and boundary conditions of a
    Nek5000/nekRS .re2 file written by gmsh2nek ('MSH' faces carry the
    physical-group ID in bc(5), 'P  ' faces the neighbour element/face in
    bc(1)/bc(2))."""
    raw = open(path, "rb").read()
    hdr = raw[:80].decode(errors="replace").split()
    nel, ndim = int(hdr[1]), int(hdr[2])
    if ndim != 3:
        raise ValueError(f"{path}: expected a 3D mesh, got ndim = {ndim}")
    bo = "<" if abs(np.frombuffer(raw[80:84], "<f4")[0] - 6.54321) < 1e-4 else ">"
    dbl = np.dtype(bo + "f8")
    off = 84
    xyz = np.frombuffer(raw, dbl, count=nel * 25, offset=off).reshape(nel, 25)
    off += nel * 25 * 8
    # Nek vertex order 1 2 4 3 / 5 6 8 7 -> tensor order i + 2j + 4k.
    sym = [0, 1, 3, 2, 4, 5, 7, 6]
    V = np.stack([xyz[:, 1:9][:, sym], xyz[:, 9:17][:, sym], xyz[:, 17:25][:, sym]], axis=2)
    dets = []
    for c in range(8):
        i, j, k = c & 1, (c >> 1) & 1, (c >> 2) & 1
        s = (-1) ** (i + j + k)
        dets.append(s * np.einsum("ij,ij->i", V[:, c ^ 1] - V[:, c], np.cross(V[:, c ^ 2] - V[:, c], V[:, c ^ 4] - V[:, c])))
    dets = np.array(dets)
    out = {"file": path, "nel": nel, "corner_jacobian_min": float(dets.min()),
           "n_nonpositive_corner_jacobian": int(np.sum(dets.min(axis=0) <= 0))}
    ncurve = int(np.frombuffer(raw, dbl, count=1, offset=off)[0])
    off += 8
    cur = np.frombuffer(raw, dbl, count=ncurve * 8, offset=off).reshape(ncurve, 8)
    off += ncurve * 64
    ctypes = [bytes(r[7:8].view("S8")[0])[:1].decode() for r in cur]
    out["curved_edges"] = {k: ctypes.count(k) for k in sorted(set(ctypes))}
    out["n_elements_with_curved_edges"] = int(len(np.unique(cur[:, 0]))) if ncurve else 0
    nbc = int(np.frombuffer(raw, dbl, count=1, offset=off)[0])
    off += 8
    b = np.frombuffer(raw, dbl, count=nbc * 8, offset=off).reshape(nbc, 8)
    off += nbc * 64
    cb = np.array([bytes(r[7:8].view("S8")[0])[:3].decode() for r in b])
    ids, n = np.unique(b[cb == "MSH", 6].astype(int), return_counts=True)
    out["msh_faces_by_id"] = {int(i): int(k) for i, k in zip(ids, n)}
    out["n_periodic_faces"] = int(np.sum(cb == "P  "))
    out["other_bc_types"] = sorted(set(cb) - {"MSH", "P  "})
    links = {(int(r[0]), int(r[1]), int(r[2]), int(r[3])) for r in b[cb == "P  "]}
    out["periodic_links_asymmetric"] = sum((e2, f2, e, f) not in links for e, f, e2, f2 in links)
    out["bytes_unread"] = len(raw) - off  # further fields (none from gmsh2nek)
    return out


def evaluate_re2(r, m):
    """Target name -> (value, limit, ok) for the .re2 summary r of the mesh with metrics m."""
    g = m["boundary_groups"]
    ids = {1: g.get("inflow", 0), 2: g.get("outflow", 0)}
    n_per = sum(g.get(k, 0) for k in ("xmin", "xmax", "zmin", "zmax"))
    return {"re2_elements": (r["nel"], m["n_elements"], r["nel"] == m["n_elements"]),
            "re2_corner_jacobians": (r["n_nonpositive_corner_jacobian"], 0, r["n_nonpositive_corner_jacobian"] == 0),
            "re2_boundary_ids": (r["msh_faces_by_id"], ids, r["msh_faces_by_id"] == ids and not r["other_bc_types"]),
            "re2_periodic_faces": (r["n_periodic_faces"], n_per, r["n_periodic_faces"] == n_per),
            "re2_periodic_links_symmetric": (r["periodic_links_asymmetric"], 0, r["periodic_links_asymmetric"] == 0),
            "re2_fully_read": (r["bytes_unread"], 0, r["bytes_unread"] == 0)}


def table(res):
    """Printable table of evaluate() results."""
    def fmt(v):
        return f"{v:.4g}" if isinstance(v, float) else str(v)
    return "\n".join(f"  {'ok  ' if ok else 'FAIL'} {name:30s} {fmt(v):>12s}  (limit {fmt(lim)})"
                     for name, (v, lim, ok) in res.items())


def summary(m):
    """A few lines on size, quality and the nekRS-LS settings of the metrics m."""
    q, iz = m["quality"], m["regions"]["interface_zone"]
    lines = [f"elements {m['n_elements']}, GLL points at N = 7: {m['gll_points_N7']}",
             f"minSJ {q['minSJ']['min']:.3f} at {q['worst_minSJ'][0]['centroid']}; corner SJ "
             f"{q['corner_SJ']['min']:.3f} at {q['corner_SJ']['worst_at']['centroid']}; CFL length "
             f"{q['cfl_length']['min']:.4f} at {q['cfl_length']['worst_at']['centroid']}",
             f"neighbour h_J ratio {m['neighbor_hJ_ratio']['max']:.3f}, thickness ratio "
             f"{m['neighbor_thickness_ratio']['max']:.3f} at {m['neighbor_thickness_ratio']['worst_at']}; "
             f"edge aspect max {q['edge_aspect']['max']:.2f}"]
    if iz:
        lines.append(f"interface zone: {iz['n_elem']} elements, edges <= {iz['edge_max']['value']:.4f}, "
                     f"h_J {iz['hJ_min']:.4f} .. {iz['hJ_max']:.4f} (median {iz['hJ_median']:.4f})")
    for key, s in m.get("nekrs_ls", {}).items():
        lines.append(f"nekRS-LS {key}: interfaceWidthValue {s['recommended_interfaceWidthValue']:.4f} "
                     f"(R/eps {s['R_over_eps']:.1f}), capillary dt {s['capillary_dt_0.8x']:.3g}, advective dt "
                     f"{s['advective_dt_cfl0.5_U1']:.3g} (CFL 0.5, U = 1), Hmax/Hmin {s['Hmax_over_Hmin']:.2f}")
    return "\n".join(lines)


def limit_memory(gb):
    """Cap this process's address space at `gb` GiB (0: no cap; a lower
    existing limit stays in force), so a runaway mesh fails instead of
    exhausting the machine (on WSL an out-of-memory kill of any process stops
    the whole distro). The cap counts virtual address space, about 2.2 GB of
    which the gmsh and numpy imports take; when it is hit, numpy raises
    MemoryError, but native code (gmsh) may abort the process instead."""
    if not 0 <= gb < math.inf:
        raise ValueError(f"max_memory_gb = {gb} must be finite and >= 0 (0: no cap)")
    if gb > 0:
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        lim = min([int(gb * 2 ** 30)] + [v for v in (soft, hard) if v != resource.RLIM_INFINITY])
        resource.setrlimit(resource.RLIMIT_AS, (lim, hard))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("msh")
    ap.add_argument("--re2", help="gmsh2nek output to check as well")
    ap.add_argument("--json", help="write the metrics and results here")
    ap.add_argument("--plan", help="read the targets from generate_mesh.py's NAME.plan.json")
    ap.add_argument("--max_memory_gb", type=float, default=8.0, help="address-space cap (0: none), default 8")
    for k, v in TARGETS.items():
        ap.add_argument(f"--{k}", type=type(v))
    a = ap.parse_args()
    limit_memory(a.max_memory_gb)
    t = dict(TARGETS)
    if a.plan:
        with open(a.plan) as f:
            t.update({k: v for k, v in json.load(f)["parameters"].items() if k in TARGETS})
    t.update({k: getattr(a, k) for k in TARGETS if getattr(a, k) is not None})
    m = measure(load_msh(a.msh), t)
    res = evaluate(m, t)
    if a.re2:
        m["re2"] = read_re2(a.re2)
        res.update(evaluate_re2(m["re2"], m))
    m["targets"] = t
    m["results"] = {k: {"value": v, "limit": lim, "ok": ok} for k, (v, lim, ok) in res.items()}
    if a.json:
        with open(a.json, "w") as f:
            json.dump(m, f, indent=1, default=str)
    print(f"{a.msh}{' + ' + a.re2 if a.re2 else ''}:\n{summary(m)}\n{table(res)}")
    failed = [k for k, (_, _, ok) in res.items() if not ok]
    if failed:
        print("FAILED:", ", ".join(failed))
        sys.exit(1)


if __name__ == "__main__":
    main()
