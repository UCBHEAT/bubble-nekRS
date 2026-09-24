"""Postprocess and animate one run of this study with ParaView.

  pvbatch animate.py post <run dir> <first> <last>
      Bubble statistics for checkpoints first..last-1 (0-based) ->
      <run dir>/post/post_<first>.csv
  python3 animate.py merge <run dir>
      Join the chunks, unwrap the periodic centroid track -> <run dir>/post.csv
  pvbatch animate.py frames <run dir> <first> <last> <label>
      Render checkpoints first..last-1 -> <run dir>/frames/frame_<i>.png

The look follows answinter26/example/paraview.pvsm: domain outline, a z = 0
slice colored by c on a log scale, the psi = 0.5 bubble surface, velocity
arrows on the slice, and bubble rise velocity, MTC and Sh charts with a
current time marker. pvbatch needs the EGL ParaView on a GPU node.
"""
import builtins
import csv
import glob
import math
import os
import sys

import numpy as np

# 4x8x4 fully periodic domain (bubble3d.box).
LO = np.array([-2.0, 0.0, -2.0])
L = np.array([4.0, 8.0, 4.0])
T_END = 30.0

POST_COLUMNS = [
    "Time", "gas_volume", "d_eq", "x_c", "y_c", "z_c",
    "u_b", "v_b", "w_b", "u_l", "v_l", "w_l", "rise_velocity", "slip_speed",
    "contour_area", "width_x", "width_z", "height_y", "aspect_ratio",
    "max_speed_gas", "max_speed_liquid", "c_gas_mean",
]


def checkpoints(run):
    """Every checkpoint of a run as (part .nek5000 file, time), in time order.

    Each job's checkpoints (part<N>/, see prepare-restart.sh) are written in
    that job's element order with the coordinates in its first file, so every
    part is read as its own series rather than through the run's links.
    """
    parts = glob.glob(os.path.join(run, "part[0-9]*"))
    parts = sorted(parts, key=lambda d: int(os.path.basename(d)[4:])) or [run]
    out = []
    for part in parts:
        files = sorted(glob.glob(os.path.join(part, "bubble3d0.f[0-9][0-9][0-9][0-9][0-9]")))
        if not files:
            continue
        nek = os.path.join(part, "animate.nek5000")
        tmp = "{}.{}".format(nek, os.getpid())
        with open(tmp, "w") as f:
            f.write(" filetemplate: bubble3d%01d.f%05d\n firsttimestep: 0\n"
                    " numtimesteps: {}\n".format(len(files)))
        os.rename(tmp, nek)
        for fn in files:
            with open(fn, "rb") as fh:
                out.append((nek, float(fh.read(132).split()[7])))
    return sorted(out, key=lambda c: c[1])


def open_reader(nek):
    from paraview.simple import Nek5000Reader
    r = Nek5000Reader(FileName=nek)
    r.PointArrays = ["Velocity", "S02", "S03"]  # S01 = tls, S02 = cls (psi), S03 = c
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


def circular_mean(weights, cos_theta, sin_theta, lo, length):
    theta = math.atan2(weights @ sin_theta, weights @ cos_theta) % (2*math.pi)
    return lo + length*theta/(2*math.pi)


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
            # The mesh is static: weights and trig tables are reused every step.
            W = point_weights(grid)
            xyz = vtk_to_numpy(grid.GetPoints().GetData()).astype(np.float64)
            theta = 2*np.pi*(xyz - LO)/L
            cos_t, sin_t = np.cos(theta), np.sin(theta)
            del theta
        pd = grid.GetPointData()
        U = vtk_to_numpy(pd.GetArray("Velocity")).astype(np.float64)
        psi = np.clip(vtk_to_numpy(pd.GetArray("S02")).astype(np.float64), 0.0, 1.0)
        c = vtk_to_numpy(pd.GetArray("S03")).astype(np.float64)
        wg = W*(1.0 - psi)
        wl = W*psi
        vg, vl = wg.sum(), wl.sum()
        center = np.array([circular_mean(wg, cos_t[:, k], sin_t[:, k], LO[k], L[k]) for k in range(3)])
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
            # Extents relative to the centroid, across periodic boundaries.
            rel = (p - center + L/2) % L - L/2
            ext = rel.max(axis=0) - rel.min(axis=0)
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
        w = csv.DictWriter(f, fieldnames=POST_COLUMNS)
        w.writeheader()
        w.writerows(rows)
    os.rename(out + ".tmp", out)


def merge(run):
    rows = []
    for fn in sorted(glob.glob(os.path.join(run, "post", "post_*.csv"))):
        with open(fn) as f:
            rows += [{k: float(v) for k, v in row.items()} for row in csv.DictReader(f)]
    rows.sort(key=lambda row: row["Time"])
    # Unwrap the centroid track so the bubble's distance travelled is continuous.
    prev = None
    for row in rows:
        now = np.array([row["x_c"], row["y_c"], row["z_c"]])
        if prev is None:
            unwrapped = now.copy()
        else:
            unwrapped = unwrapped + (now - prev + L/2) % L - L/2
        prev = now
        row["x_c_unwrapped"], row["y_c_unwrapped"], row["z_c_unwrapped"] = unwrapped
    columns = POST_COLUMNS + ["x_c_unwrapped", "y_c_unwrapped", "z_c_unwrapped"]
    with open(os.path.join(run, "post.csv"), "w") as f:
        w = csv.DictWriter(f, fieldnames=columns)
        w.writeheader()
        for row in rows:
            w.writerow({k: "{:.6g}".format(row[k]) for k in columns})
    print("{}: merged {} rows into post.csv".format(run, len(rows)))


def frames(run, first, last, label):
    from paraview.simple import (CSVReader, ColorBy, Contour, CreateLayout, CreateView, Delete,
                                 GetColorTransferFunction, GetScalarBar, Glyph, Hide,
                                 ProgrammableSource, SaveScreenshot, Show, Slice, Text)
    cps = checkpoints(run)
    current = cps[builtins.min(first, len(cps) - 1)][0]
    r = open_reader(current)
    os.makedirs(os.path.join(run, "frames"), exist_ok=True)

    rv = CreateView("RenderView")
    rv.OrientationAxesVisibility = 1
    rv.UseColorPaletteForBackground = 0
    rv.Background = [0.32, 0.34, 0.43]
    # Camera from the example state, scaled from its 2x4x2 domain to 4x8x4.
    focal = np.array([0.0, 4.0, 0.0])
    rv.CameraFocalPoint = list(focal)
    rv.CameraPosition = list(focal + 2.0*np.array([5.795554957734411, 4.7320508075688785, 5.795554957734411]))
    rv.CameraViewUp = [-0.35355339059327373, 0.8660254037844388, -0.35355339059327373]
    rv.CameraViewAngle = 30

    def show_outline(src):
        disp = Show(src, rv)
        disp.SetRepresentationType("Outline")
        disp.ColorArrayName = ["POINTS", ""]
        disp.AmbientColor = disp.DiffuseColor = [1, 1, 1]

    show_outline(r)

    slc = Slice(Input=r)
    slc.SliceType = "Plane"
    slc.SliceType.Origin = list(focal)
    slc.SliceType.Normal = [0, 0, 1]
    slc_disp = Show(slc, rv)
    ColorBy(slc_disp, ("POINTS", "S03"))
    c_lut = GetColorTransferFunction("S03")
    c_lut.AutomaticRescaleRangeMode = "Never"
    c_lut.UseLogScale = 1
    c_lut.RGBPoints = [1e-4, 1, 1, 1, 1.0, 0, 0.3333333333333333, 1]
    c_lut.RescaleTransferFunction(1e-4, 1.0)
    slc_disp.SetScalarBarVisibility(rv, True)
    c_bar = GetScalarBar(c_lut, rv)
    c_bar.Title = "c"
    c_bar.ComponentTitle = ""
    c_bar.WindowLocation = "Lower Right Corner"

    bubble = Contour(Input=r, ContourBy=["POINTS", "S02"], Isosurfaces=[0.5])
    bubble_disp = Show(bubble, rv)
    bubble_disp.ColorArrayName = ["POINTS", ""]
    bubble_disp.AmbientColor = bubble_disp.DiffuseColor = [0, 1, 1]
    bubble_disp.Opacity = 0.5

    glyph = Glyph(Input=slc, GlyphType="Arrow")
    glyph.OrientationArray = ["POINTS", "Velocity"]
    glyph.ScaleArray = ["POINTS", "Velocity"]
    glyph.ScaleFactor = 0.4
    glyph.GlyphMode = "Uniform Spatial Distribution (Surface Sampling)"
    glyph.MaximumNumberOfSamplePoints = 5000
    glyph_disp = Show(glyph, rv)
    ColorBy(glyph_disp, ("POINTS", "Velocity", "Magnitude"))
    u_lut = GetColorTransferFunction("Velocity")
    u_lut.AutomaticRescaleRangeMode = "Never"
    u_lut.RGBPoints = [0, 1, 1, 1, 3, 1, 0, 0]
    u_lut.RescaleTransferFunction(0.0, 3.0)
    glyph_disp.SetScalarBarVisibility(rv, True)
    u_bar = GetScalarBar(u_lut, rv)
    u_bar.Title = "velocity magnitude"
    u_bar.ComponentTitle = ""
    u_bar.WindowLocation = "Upper Right Corner"
    for bar in (c_bar, u_bar):
        bar.TitleFontSize = 24
        bar.LabelFontSize = 22
        bar.ScalarBarLength = 0.3

    label_text = Text()
    label_disp = Show(label_text, rv)
    label_disp.WindowLocation = "Upper Left Corner"
    label_disp.FontSize = 26

    data_csv = CSVReader(FileName=[os.path.join(run, "data.csv")])
    post_csv = CSVReader(FileName=[os.path.join(run, "post.csv")])
    ymax = {"vel": 1.5, "mtc": 0.2, "sh": 40.0}
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

    def chart(title, src, series, color, ymax, marker_series, bottom_title=""):
        ch = CreateView("XYChartView")
        ch.ChartTitle = title
        ch.ChartTitleFontSize = 30
        ch.ShowLegend = 0
        ch.LeftAxisUseCustomRange = 1
        ch.LeftAxisRangeMinimum, ch.LeftAxisRangeMaximum = 0.0, ymax
        ch.BottomAxisUseCustomRange = 1
        ch.BottomAxisRangeMinimum, ch.BottomAxisRangeMaximum = 0.0, T_END
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
    vel_chart = chart("Bubble rise velocity", post_csv, "rise_velocity", green, ymax["vel"], "vel")
    mtc_chart = chart("Mass transfer coefficient", data_csv, "MTC", red, ymax["mtc"], "mtc")
    sh_chart = chart("Sherwood number", data_csv, "Sh", green, ymax["sh"], "sh", "Time")

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
            old, r = r, open_reader(nek)
            slc.Input = r
            bubble.Input = r
            Hide(old, rv)
            show_outline(r)
            Delete(old)
            current = nek
        set_marker(t)
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
