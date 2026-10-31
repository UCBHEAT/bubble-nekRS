"""Postprocess and animate one run of this study with ParaView.

gpu/mesh-sensitivity's animate.py for the PID-centred case: the run's file
names come from its .par, the bubble stays near the middle of the (only x
and z periodic) domain, so plain averages and extents replace the periodic
ones, the view follows the bubble and its wake, the arrows show the velocity
relative to the far-field liquid (the frame velocity from data.csv
subtracted), c is on a linear scale (it stays near 1 away from the wake)
and the rise velocity chart is data.csv's lab-frame rise_v.

  pvbatch animate.py post <run dir> <first> <last>
      Bubble statistics for checkpoints first..last-1 (0-based) ->
      <run dir>/post/post_<first>.csv
  python3 animate.py merge <run dir>
      Join the chunks, unwrap the periodic centroid track -> <run dir>/post.csv
      (with the liquid statistics of <run dir>/full-domain.csv, if any: see
      subset-checkpoints.py)
  pvbatch animate.py frames <run dir> <first> <last> <label>
      Render checkpoints first..last-1 -> <run dir>/frames/frame_<i>.png

The look follows answinter26/example/paraview.pvsm: a z = 0 slice colored
by c, the psi = 0.5 bubble surface, velocity arrows on the slice, and bubble
rise velocity, MTC and Sh charts with a current time marker. pvbatch needs
the EGL ParaView on a GPU node.

ANIMATE_SPECIES=c<Sc> shows a further species of the run instead of c: its
field on the slice and its MTC and Sh (data-c<Sc>.csv) in the charts. The
charts' time axes start where the run started (data.csv's first row): at 0,
or at START for a run that setup-runs.sh started FROM another.
"""
import builtins
import csv
import glob
import math
import os
import sys

import numpy as np

T_END = 30.0  # chart time axis when the run's .par has no endTime
SPECIES = os.environ.get("ANIMATE_SPECIES", "c")

POST_COLUMNS = [
    "Time", "gas_volume", "d_eq", "x_c", "y_c", "z_c",
    "u_b", "v_b", "w_b", "u_l", "v_l", "w_l", "rise_velocity", "slip_speed",
    "contour_area", "width_x", "width_z", "height_y", "aspect_ratio",
    "max_speed_gas", "max_speed_liquid", "c_gas_mean",
]


def case_name(run):
    """The run's case name: that of its only .par file (bubble here)."""
    pars = glob.glob(os.path.join(run, "*.par"))
    assert len(pars) == 1, "{} must hold exactly one .par file".format(run)
    return os.path.basename(pars[0])[:-len(".par")]


def checkpoints(run):
    """Every checkpoint of a run as (part .nek5000 file, time), in time order.

    Each job's checkpoints (part<N>/, see prepare-restart.sh) are written in
    that job's element order with the coordinates in its first file, so every
    part is read as its own series rather than through the run's links.
    """
    parts = glob.glob(os.path.join(run, "part[0-9]*"))
    parts = sorted(parts, key=lambda d: int(os.path.basename(d)[4:])) or [run]
    case = case_name(run)
    out = []
    for part in parts:
        files = sorted(glob.glob(os.path.join(part, case + "0.f[0-9][0-9][0-9][0-9][0-9]")))
        if not files:
            continue
        nek = os.path.join(part, "animate.nek5000")
        tmp = "{}.{}".format(nek, os.getpid())
        with open(tmp, "w") as f:
            f.write(" filetemplate: {}%01d.f%05d\n firsttimestep: 0\n"
                    " numtimesteps: {}\n".format(case, len(files)))
        os.rename(tmp, nek)
        for fn in files:
            with open(fn, "rb") as fh:
                out.append((nek, float(fh.read(132).split()[7])))
    return sorted(out, key=lambda c: c[1])


def end_time(run):
    """The run's endTime from its .par (T_END if it has none), the end of the
    chart time axes."""
    with open(os.path.join(run, case_name(run) + ".par")) as f:
        for line in f:
            key, _, value = line.partition("=")
            if key.strip() == "endTime":
                return float(value.split("#")[0])
    return T_END


def start_time(run, data="data.csv"):
    """When the run (or the species of data) started: the time of the data
    file's first row, rounded down; 0 from rest, START for a run setup-runs.sh
    started FROM another."""
    with open(os.path.join(run, data)) as f:
        first = next(csv.DictReader(f), None)
    return float(math.floor(float(first["Time"]))) if first else 0.0


def species_field(run):
    """SPECIES's field in the checkpoints, S<n> for the n-th of the .par's
    scalars (S01 = tls, S02 = cls (psi), S03 = c, then the further species),
    and its data file."""
    names = ["tls", "cls", "c"]
    with open(os.path.join(run, case_name(run) + ".par")) as f:
        for line in f:
            key, _, value = line.partition("=")
            if key.strip() == "scalars":
                names = [s.strip() for s in value.split("#")[0].split(",")]
    if SPECIES not in names:
        sys.exit("{}: no scalar {} in its .par".format(run, SPECIES))
    data = "data.csv" if SPECIES == "c" else "data-{}.csv".format(SPECIES)
    return "S{:02d}".format(names.index(SPECIES) + 1), data


def chart_max(run, column, default, tmin=1.0, data="data.csv"):
    """Chart axis maximum, shared by the runs of a study (the run directories
    next to this one): the default, unless their data after the initial
    transient (the first tmin time units of each) exceed it, which it would
    clip, or stay below a fifth of it (another Schmidt number or case); then
    the round number (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6 or 8 times a power of
    10) just above 1.1 times their largest value."""
    top = None
    for path in glob.glob(os.path.join(os.path.dirname(os.path.abspath(run)), "*", data)):
        t0 = start_time(os.path.dirname(path), data) + tmin
        with open(path) as f:
            values = [float(r[column]) for r in csv.DictReader(f) if float(r["Time"]) > t0]
        if values:
            top = builtins.max(values) if top is None else builtins.max([top] + values)
    if top is None or top <= 0 or 0.2*default <= top <= default:
        return default
    scale = 10.0**math.floor(math.log10(1.1*top))
    return builtins.min(m*scale for m in (1, 1.2, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10) if m*scale >= 1.1*top)


def frame_velocities(run):
    """The frame velocity (data.csv frame_u, frame_v, frame_w) at a time: the
    row nearest to it (checkpoints fall on rows)."""
    with open(os.path.join(run, "data.csv")) as f:
        rows = [(float(r["Time"]), float(r["frame_u"]), float(r["frame_v"]), float(r["frame_w"]))
                for r in csv.DictReader(f)]
    return lambda t: builtins.min(rows, key=lambda row: abs(row[0] - t))[1:]


def open_reader(nek, field="S03"):
    from paraview.simple import Nek5000Reader
    r = Nek5000Reader(FileName=nek)
    r.PointArrays = ["Velocity", "S02", field]  # S01 = tls, S02 = cls (psi), S03 = c
    r.UpdatePipelineInformation()
    return r


def single_block(obj):
    """The reader and contour may hand back a multiblock with one leaf."""
    from vtkmodules.vtkCommonDataModel import vtkCompositeDataSet
    if isinstance(obj, vtkCompositeDataSet):
        it = obj.NewIterator()
        it.InitTraversal()
        leaves = []
        while not it.IsDoneWithTraversal():
            leaves.append(it.GetCurrentDataObject())
            it.GoToNextItem()
        assert len(leaves) == 1, "expected a single block"
        return leaves[0]
    return obj


def point_weights(grid):
    """Quadrature weight of each point: 1/8 of every adjacent hex volume.

    The box mesh hexes are axis aligned, so each volume is the product of its
    corners' extents.
    """
    from vtkmodules.util.numpy_support import vtk_to_numpy
    xyz = vtk_to_numpy(grid.GetPoints().GetData())
    conn = vtk_to_numpy(grid.GetCells().GetConnectivityArray()).reshape(-1, 8)
    vol = np.ones(len(conn))
    for k in range(3):
        corner = xyz[:, k][conn]
        vol *= corner.max(axis=1) - corner.min(axis=1)
    return np.bincount(conn.ravel(), weights=np.repeat(vol/8.0, 8), minlength=len(xyz))


def post(run, first, last):
    from paraview.simple import Contour, Delete
    from vtkmodules.util.numpy_support import vtk_to_numpy
    os.makedirs(os.path.join(run, "post"), exist_ok=True)
    out = os.path.join(run, "post", "post_{:05d}.csv".format(first))
    if os.path.exists(out):
        print("{} exists, skipping".format(out))
        return
    cps = checkpoints(run)
    r = contour = current = None
    rows = []
    for i in range(first, builtins.min(last, len(cps))):
        nek, t = cps[i]
        if nek != current:
            # A fresh reader per part: the Nek5000 reader keeps serving the old
            # part's data after a FileName change.
            if r is not None:
                Delete(contour)
                Delete(r)
            r = open_reader(nek)
            contour = Contour(Input=r, ContourBy=["POINTS", "S02"], Isosurfaces=[0.5],
                              ComputeNormals=0, ComputeScalars=0)
            current = nek
            W = None  # each part has its own element (and point) order
        r.UpdatePipeline(t)
        grid = single_block(r.GetClientSideObject().GetOutputDataObject(0))
        if W is None:
            # The mesh is static: the weights are reused every step.
            W = point_weights(grid)
            xyz = vtk_to_numpy(grid.GetPoints().GetData()).astype(np.float64)
        pd = grid.GetPointData()
        U = vtk_to_numpy(pd.GetArray("Velocity")).astype(np.float64)
        psi = np.clip(vtk_to_numpy(pd.GetArray("S02")).astype(np.float64), 0.0, 1.0)
        c = vtk_to_numpy(pd.GetArray("S03")).astype(np.float64)
        wg = W*(1.0 - psi)
        wl = W*psi
        vg, vl = wg.sum(), wl.sum()
        # The PID controller keeps the bubble near the middle of the domain,
        # far from the periodic x and z boundaries, so plain means do.
        center = wg @ xyz/vg
        ub = wg @ U/vg
        ul = wl @ U/vl
        speed = np.sqrt((U*U).sum(axis=1))

        contour.UpdatePipeline(t)
        surf = single_block(contour.GetClientSideObject().GetOutputDataObject(0))
        if surf.GetNumberOfPoints() > 0:
            p = vtk_to_numpy(surf.GetPoints().GetData()).astype(np.float64)
            tri = vtk_to_numpy(surf.GetPolys().GetConnectivityArray()).reshape(-1, 3)
            area = 0.5*np.linalg.norm(np.cross(p[tri[:, 1]] - p[tri[:, 0]],
                                               p[tri[:, 2]] - p[tri[:, 0]]), axis=1).sum()
            ext = p.max(axis=0) - p.min(axis=0)
        else:
            area, ext = 0.0, np.zeros(3)
        rows.append({
            "Time": t, "gas_volume": vg, "d_eq": (6.0*vg/math.pi)**(1.0/3.0),
            "x_c": center[0], "y_c": center[1], "z_c": center[2],
            "u_b": ub[0], "v_b": ub[1], "w_b": ub[2],
            "u_l": ul[0], "v_l": ul[1], "w_l": ul[2],
            "rise_velocity": ub[1] - ul[1], "slip_speed": np.linalg.norm(ub - ul),
            "contour_area": area, "width_x": ext[0], "width_z": ext[2], "height_y": ext[1],
            "aspect_ratio": 0.5*(ext[0] + ext[2])/ext[1] if ext[1] > 0 else float("nan"),
            "max_speed_gas": speed[psi < 0.5].max(), "max_speed_liquid": speed[psi > 0.95].max(),
            "c_gas_mean": (wg @ c)/vg,
        })
        print("post {} t={:.4f} V_g={:.4f} y_c={:.3f} rise={:.3f} area={:.3f}".format(
            os.path.basename(run), t, vg, center[1], ub[1] - ul[1], area), flush=True)
    with open(out + ".tmp", "w") as f:
        w = csv.DictWriter(f, fieldnames=POST_COLUMNS, lineterminator="\n")
        w.writeheader()
        w.writerows(rows)
    os.rename(out + ".tmp", out)


def merge(run):
    rows = []
    for fn in sorted(glob.glob(os.path.join(run, "post", "post_*.csv"))):
        with open(fn) as f:
            rows += [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
    rows.sort(key=lambda row: row["Time"])
    # A run animated from subset-checkpoints.py's cut-down checkpoints: the
    # liquid statistics span the whole domain, which the subset lacks.
    full = os.path.join(run, "full-domain.csv")
    if os.path.exists(full):
        with open(full) as f:
            whole = [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
        for row in rows:
            match = [w for w in whole if abs(w["Time"] - row["Time"]) < 1e-9]
            if len(match) != 1:
                sys.exit("{}: no single row for t = {}".format(full, row["Time"]))
            for k in ("u_l", "v_l", "w_l", "max_speed_liquid"):
                row[k] = match[0][k]
            row["rise_velocity"] = row["v_b"] - row["v_l"]
            row["slip_speed"] = math.sqrt(sum((row[b] - row[l])**2 for b, l in
                                              (("u_b", "u_l"), ("v_b", "v_l"), ("w_b", "w_l"))))
    columns = POST_COLUMNS
    with open(os.path.join(run, "post.csv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=columns, lineterminator="\n")
        w.writeheader()
        for row in rows:
            w.writerow({k: "{:.6g}".format(row[k]) for k in columns})
    print("{}: merged {} rows into post.csv".format(run, len(rows)))


def frames(run, first, last, label):
    from paraview.simple import (CSVReader, Calculator, ColorBy, Contour, CreateLayout, CreateView,
                                 Delete, GetColorTransferFunction, GetScalarBar, Glyph,
                                 ProgrammableSource, SaveScreenshot, Show, Slice, Text)
    cps = checkpoints(run)
    field, data = species_field(run)
    current = cps[builtins.min(first, len(cps) - 1)][0]
    r = open_reader(current, field)
    os.makedirs(os.path.join(run, "frames"), exist_ok=True)

    rv = CreateView("RenderView")
    rv.OrientationAxesVisibility = 1
    rv.UseColorPaletteForBackground = 0
    # White, with dark text: the saturated liquid of the slice is near white.
    rv.Background = [1.0, 1.0, 1.0]
    dark = [0.1, 0.1, 0.1]
    rv.OrientationAxesLabelColor = dark
    # Camera from the example state, scaled from its 2x4x2 domain to show the
    # bubble (at the origin) and about 5 D of its wake, not the whole 6x15x6
    # domain.
    scale = 1.5
    focal = np.array([0.0, -1.5, 0.0])
    rv.CameraFocalPoint = list(focal)
    rv.CameraPosition = list(focal + scale*np.array([5.795554957734411, 4.7320508075688785, 5.795554957734411]))
    rv.CameraViewUp = [-0.35355339059327373, 0.8660254037844388, -0.35355339059327373]
    rv.CameraViewAngle = 30

    slc = Slice(Input=r)
    slc.SliceType = "Plane"
    slc.SliceType.Origin = [0.0, 0.0, 0.0]
    slc.SliceType.Normal = [0, 0, 1]
    slc_disp = Show(slc, rv)
    # Unlit, so the oblique slice shows the colormap's own colors.
    slc_disp.Ambient, slc_disp.Diffuse = 1.0, 0.0
    ColorBy(slc_disp, ("POINTS", field))
    c_lut = GetColorTransferFunction(field)
    c_lut.AutomaticRescaleRangeMode = "Never"
    # Light where the liquid is saturated, dark in the depleted wake.
    c_lut.UseLogScale = 0
    c_lut.RGBPoints = [0.0, 0.031, 0.188, 0.42, 0.5, 0.42, 0.68, 0.84, 1.0, 0.97, 0.98, 1.0]
    c_lut.RescaleTransferFunction(0.0, 1.0)
    slc_disp.SetScalarBarVisibility(rv, True)
    c_bar = GetScalarBar(c_lut, rv)
    c_bar.Title = SPECIES
    c_bar.ComponentTitle = ""
    c_bar.WindowLocation = "Lower Right Corner"

    bubble = Contour(Input=r, ContourBy=["POINTS", "S02"], Isosurfaces=[0.5])
    bubble_disp = Show(bubble, rv)
    bubble_disp.ColorArrayName = ["POINTS", ""]
    bubble_disp.AmbientColor = bubble_disp.DiffuseColor = [0, 1, 1]
    bubble_disp.Opacity = 0.5

    # The velocity relative to the far-field liquid, which moves at the frame
    # velocity of data.csv: the flow the bubble induces, as in a lab frame.
    frame_velocity = frame_velocities(run)
    lab = Calculator(Input=slc, ResultArrayName="u_lab")

    def set_frame(t):
        fu, fv, fw = frame_velocity(t)
        lab.Function = "Velocity - ({:.10g}*iHat + {:.10g}*jHat + {:.10g}*kHat)".format(fu, fv, fw)

    set_frame(cps[builtins.min(first, len(cps) - 1)][1])
    glyph = Glyph(Input=lab, GlyphType="Arrow")
    glyph.OrientationArray = ["POINTS", "u_lab"]
    glyph.ScaleArray = ["POINTS", "u_lab"]
    glyph.ScaleFactor = 0.2*scale
    glyph.GlyphMode = "Uniform Spatial Distribution (Surface Sampling)"
    glyph.MaximumNumberOfSamplePoints = 5000
    glyph_disp = Show(glyph, rv)
    ColorBy(glyph_disp, ("POINTS", "u_lab", "Magnitude"))
    u_lut = GetColorTransferFunction("u_lab")
    u_lut.AutomaticRescaleRangeMode = "Never"
    u_lut.RGBPoints = [0, 0.75, 0.75, 0.75, 2, 1, 0, 0]
    u_lut.RescaleTransferFunction(0.0, 2.0)
    glyph_disp.SetScalarBarVisibility(rv, True)
    u_bar = GetScalarBar(u_lut, rv)
    u_bar.Title = "velocity relative to the far field"
    u_bar.ComponentTitle = ""
    u_bar.WindowLocation = "Upper Right Corner"
    for bar in (c_bar, u_bar):
        bar.TitleFontSize = 24
        bar.LabelFontSize = 22
        bar.TitleColor = bar.LabelColor = dark
        bar.ScalarBarLength = 0.3

    label_text = Text()
    label_disp = Show(label_text, rv)
    label_disp.WindowLocation = "Upper Left Corner"
    label_disp.FontSize = 26
    label_disp.Color = dark

    data_csv = CSVReader(FileName=[os.path.join(run, "data.csv")])
    species_csv = data_csv if data == "data.csv" else CSVReader(FileName=[os.path.join(run, data)])
    ymax = {"vel": chart_max(run, "rise_v", 1.5),
            "mtc": chart_max(run, "MTC", 0.2, data=data),
            "sh": chart_max(run, "Sh", 40.0, data=data)}
    marker = ProgrammableSource(OutputDataSetType="vtkTable")

    def set_marker(t):
        marker.Script = (
            "import numpy as np\n"
            "output.RowData.append(np.array([{t}, {t}]), 'Time')\n"
            "output.RowData.append(np.array([0.0, {v}]), 'vel')\n"
            "output.RowData.append(np.array([0.0, {m}]), 'mtc')\n"
            "output.RowData.append(np.array([0.0, {s}]), 'sh')\n"
        ).format(t=t, v=ymax["vel"], m=ymax["mtc"], s=ymax["sh"])

    set_marker(0.0)
    marker.UpdatePipeline()

    t_start = start_time(run)
    # The last checkpoint can lie a few steps past endTime (13.000003 on 0.132x).
    t_end = builtins.max(end_time(run), math.ceil(cps[-1][1] - 1e-3))

    def chart(title, src, series, color, ymax, marker_series, bottom_title=""):
        ch = CreateView("XYChartView")
        ch.ChartTitle = title
        ch.ChartTitleFontSize = 30
        ch.ShowLegend = 0
        ch.LeftAxisUseCustomRange = 1
        ch.LeftAxisRangeMinimum, ch.LeftAxisRangeMaximum = 0.0, ymax
        ch.BottomAxisUseCustomRange = 1
        ch.BottomAxisRangeMinimum, ch.BottomAxisRangeMaximum = t_start, t_end
        ch.BottomAxisTitle = bottom_title
        ch.LeftAxisLabelFontSize = ch.BottomAxisLabelFontSize = 20
        ch.BottomAxisTitleFontSize = 22
        d = Show(src, ch, "XYChartRepresentation")
        d.UseIndexForXAxis = 0
        d.XArrayName = "Time"
        d.SeriesVisibility = [series]
        d.SeriesColor = [series] + [str(x) for x in color]
        d.SeriesLineThickness = [series, "2"]
        m = Show(marker, ch, "XYChartRepresentation")
        m.UseIndexForXAxis = 0
        m.XArrayName = "Time"
        m.SeriesVisibility = [marker_series]
        m.SeriesColor = [marker_series, "0.4", "0.4", "0.4"]
        m.SeriesLineThickness = [marker_series, "2"]
        return ch

    green, red = (0.3, 0.69, 0.29), (0.89, 0.1, 0.11)
    vel_chart = chart("Bubble rise velocity", data_csv, "rise_v", green, ymax["vel"], "vel")
    mtc_chart = chart("Mass transfer coefficient", species_csv, "MTC", red, ymax["mtc"], "mtc")
    sh_chart = chart("Sherwood number", species_csv, "Sh", green, ymax["sh"], "sh", "Time")

    # Render view on the left half; the three charts stacked 1:1:2 on the right.
    layout = CreateLayout("animation")
    layout.SplitHorizontal(0, 0.5)
    layout.AssignView(1, rv)
    layout.SplitVertical(2, 0.25)
    layout.AssignView(5, vel_chart)
    layout.SplitVertical(6, 1.0/3.0)
    layout.AssignView(13, mtc_chart)
    layout.AssignView(14, sh_chart)
    layout.SetSize(1618, 1000)

    # The programmable source's star import of VTK's numpy algorithms shadows
    # min/max in this namespace, hence builtins.min.
    for i in range(first, builtins.min(last, len(cps))):
        nek, t = cps[i]
        png = os.path.join(run, "frames", "frame_{:05d}.png".format(i))
        if os.path.exists(png):
            continue
        if nek != current:
            # A fresh reader per part: the Nek5000 reader keeps serving the old
            # part's data after a FileName change.
            old, r = r, open_reader(nek, field)
            slc.Input = r
            bubble.Input = r
            Delete(old)
            current = nek
        set_marker(t)
        set_frame(t)
        label_text.Text = "{}   t = {:.1f}".format(label, t)
        rv.ViewTime = t
        SaveScreenshot(png + ".tmp.png", layout, ImageResolution=[1618, 1000])
        os.rename(png + ".tmp.png", png)
        print("frame {} {} t={:.4f}".format(os.path.basename(run), i, t), flush=True)


if __name__ == "__main__":
    mode, run = sys.argv[1], sys.argv[2]
    if mode == "post":
        post(run, int(sys.argv[3]), int(sys.argv[4]))
    elif mode == "merge":
        merge(run)
    elif mode == "frames":
        frames(run, int(sys.argv[3]), int(sys.argv[4]), sys.argv[5])
    else:
        sys.exit(__doc__)
