"""ParaView animation of the PID-centered rising bubble (2D or 3D cases).

Left: 3D render of the bubble interface (cls = 0.5 contour) over the z mid-plane
slice (through the bubble centre) colored by the liquid z-vorticity, with
velocity glyphs showing the liquid flowing down past the bubble in the moving
(PID-centered) frame, and the c = 0.9 isoline of species-depleted liquid
(drawn only where psi > 0.95). Right: the PID centering force (mean over the
last 0.5 time units) and the bubble centroid offset from data.csv, drawn up to
the current animation time.

Shared by gpu/pid-centering (quasi-2D, one z element) and gpu/pid-centering-3D.
Usage, from a case directory after running it:

    ~/.local/paraview-5.13.2/bin/pvbatch ../common/animate.py [CASE_DIR] [--no-movie]

CASE_DIR defaults to the current directory (or $BUBBLE_CASE if set) and must
contain bubble.nek5000 (written by nekRS) and data.csv. Under pvbatch, frames
are written to CASE_DIR/frames/, encoded to CASE_DIR/bubble-pid.mp4 with
ffmpeg, and the pipeline is saved as CASE_DIR/animate.pvsm (open it in the
ParaView GUI via File > Load State; "Search files under specified directory"
can point it at a continued run or another case with the same domain. The
chart time axes follow the loaded data, but the camera and the velocity-glyph
grid keep the domain bounds the state was saved with, so for a different
domain rerun `pvbatch ../common/animate.py OTHER_CASE_DIR`).

In the GUI you can also run this file from View > Python Shell > Run Script to
build the same pipeline and edit it interactively (set BUBBLE_CASE, or start
ParaView from the case directory); nothing is written then.
"""

import itertools
import os
import shutil
import subprocess
import sys

from paraview.simple import *
from paraview import servermanager

args = [a for a in sys.argv[1:] if not a.startswith("--")]
CASE = os.path.abspath(args[0] if args else os.environ.get("BUBBLE_CASE", os.getcwd()))
# Only write the state, frames, and movie under pvbatch (not from the GUI).
_pm = servermanager.vtkProcessModule
IN_BATCH = _pm.GetProcessType() == _pm.PROCESS_BATCH
MAKE_MOVIE = IN_BATCH and "--no-movie" not in sys.argv

NEK = os.path.join(CASE, "bubble.nek5000")
CSV = os.path.join(CASE, "data.csv")
FRAMES = os.path.join(CASE, "frames")
MOVIE = os.path.join(CASE, "bubble-pid.mp4")
STATE = os.path.join(CASE, "animate.pvsm")

WIDTH, HEIGHT = 1920, 1080
FPS = 20
# Colormap range, fixed so colors mean the same thing in every frame. (The
# chart ranges are set from the data, see series_range.)
VORTICITY_RANGE = 5.0

# -----------------------------------------------------------------------------
# 3D data: nekRS field files. S01 = tls, S02 = cls (psi, 1 = liquid), S03 = c.
# -----------------------------------------------------------------------------
nek = Nek5000Reader(registrationName="bubble.nek5000", FileName=NEK)
nek.PointArrays = ["Velocity", "Velocity Magnitude", "Pressure", "S01", "S02", "S03"]
# Merge coincident element-boundary points so gradients are smooth across elements.
nek.MergePointstoCleangrid = 1
nek.UpdatePipelineInformation()
times = list(nek.TimestepValues) if nek.TimestepValues else [0.0]
end_time = times[-1]

nek.UpdatePipeline(times[0])
bounds = nek.GetDataInformation().GetBounds()
zmid = 0.5 * (bounds[4] + bounds[5])

vort = Gradient(registrationName="vorticity", Input=nek)
vort.ScalarArray = ["POINTS", "Velocity"]
vort.ComputeGradient = 0
vort.ComputeVorticity = 1
vort.VorticityArrayName = "Vorticity"

# Weight by the liquid fraction psi (S02) to show the liquid flow structure
# (boundary layers and wake) without the noisy gas-phase interior.
liquid_vort = Calculator(registrationName="liquid_vorticity", Input=vort)
liquid_vort.ResultArrayName = "omega_z_liquid"
liquid_vort.Function = "Vorticity_Z*S02"

midplane = Slice(registrationName="midplane", Input=liquid_vort)
midplane.SliceType = "Plane"
midplane.SliceType.Origin = [0.0, 0.0, zmid]
midplane.SliceType.Normal = [0.0, 0.0, 1.0]

# Species-depleted liquid: the c = 0.9 isoline in the mid-plane, drawn only in
# the liquid (psi > 0.95). c = 1 in the inflowing liquid and 0 in the gas (hard
# sink at psi < 0.5). Inside the interface band c rings next to that cut, so an
# unmasked isoline there traces numerical oscillations, not depletion. c starts
# as the smeared c = psi, so early on the wake is mostly that initially depleted
# band swept to the rear (and trapped in the recirculating wake), not only
# species transferred at the interface.
liquid_c = Calculator(registrationName="liquid_c", Input=midplane)
liquid_c.ResultArrayName = "c_liquid"
liquid_c.Function = "S02 > 0.95 ? S03 : 1.0"
species_wake = Contour(registrationName="species_wake", Input=liquid_c)
species_wake.ContourBy = ["POINTS", "c_liquid"]
species_wake.Isosurfaces = [0.9]

interface = Contour(registrationName="interface", Input=nek)
interface.ContourBy = ["POINTS", "S02"]
interface.Isosurfaces = [0.5]
interface.ComputeNormals = 1

# Velocity glyphs on a coarse grid of points in the mid-plane.
grid = Plane(registrationName="glyph_grid")
grid.Origin = [bounds[0], bounds[2], zmid]
grid.Point1 = [bounds[1], bounds[2], zmid]
grid.Point2 = [bounds[0], bounds[3], zmid]
# About 0.4 D between arrows.
grid.XResolution = max(8, int(round((bounds[1] - bounds[0])/0.4)))
grid.YResolution = max(8, int(round((bounds[3] - bounds[2])/0.4)))
probe = ResampleWithDataset(registrationName="glyph_probe", SourceDataArrays=nek, DestinationMesh=grid)
probe.CellLocator = "Static Cell Locator"
arrows = Glyph(registrationName="velocity_glyphs", Input=probe, GlyphType="Arrow")
arrows.OrientationArray = ["POINTS", "Velocity"]
arrows.ScaleArray = ["POINTS", "Velocity Magnitude"]
# Arrow length = ScaleFactor*|u|: ~0.35 D at |u| = 1, about the arrow spacing.
arrows.ScaleFactor = 0.35
arrows.GlyphMode = "All Points"
# Thicker than the default arrow so the shafts stay visible at this distance.
arrows.GlyphType.ShaftRadius = 0.06
arrows.GlyphType.TipRadius = 0.15

rv = CreateView("RenderView")
rv.Background = [1.0, 1.0, 1.0]
rv.UseColorPaletteForBackground = 0
rv.OrientationAxesVisibility = 0
rv.EnableRenderOnInteraction = 0

mid_disp = Show(midplane, rv)
ColorBy(mid_disp, ("POINTS", "omega_z_liquid"))
vlut = GetColorTransferFunction("omega_z_liquid")
vlut.ApplyPreset("Cool to Warm", True)
vlut.RescaleTransferFunction(-VORTICITY_RANGE, VORTICITY_RANGE)
vlut.AutomaticRescaleRangeMode = "Never"
# Unlit, so the colormap reads the same regardless of the viewing angle.
mid_disp.Ambient = 1.0
mid_disp.Diffuse = 0.0
mid_disp.SetScalarBarVisibility(rv, True)
vbar = GetScalarBar(vlut, rv)
vbar.Title = "liquid z-vorticity"
vbar.ComponentTitle = ""
vbar.TitleColor = [0.1, 0.1, 0.1]
vbar.LabelColor = [0.1, 0.1, 0.1]
vbar.TitleFontSize = 20
vbar.LabelFontSize = 16
vbar.Orientation = "Horizontal"
vbar.WindowLocation = "Lower Center"
vbar.ScalarBarLength = 0.4
vbar.RangeLabelFormat = "%-#.1f"

if_disp = Show(interface, rv)
if_disp.ColorArrayName = ["POINTS", ""]
if_disp.AmbientColor = [0.85, 0.9, 1.0]
if_disp.DiffuseColor = [0.85, 0.9, 1.0]
if_disp.Specular = 0.6
if_disp.SpecularPower = 40.0

arrow_disp = Show(arrows, rv)
arrow_disp.ColorArrayName = ["POINTS", ""]
arrow_disp.AmbientColor = [0.15, 0.15, 0.15]
arrow_disp.DiffuseColor = [0.15, 0.15, 0.15]

wake_disp = Show(species_wake, rv)
wake_disp.ColorArrayName = ["POINTS", ""]
wake_disp.AmbientColor = [0.0, 0.55, 0.25]
wake_disp.DiffuseColor = [0.0, 0.55, 0.25]
wake_disp.LineWidth = 3.0

outline_disp = Show(Outline(registrationName="domain_outline", Input=nek), rv)
outline_disp.AmbientColor = [0.3, 0.3, 0.3]
outline_disp.DiffuseColor = [0.3, 0.3, 0.3]

tlabel = AnnotateTimeFilter(registrationName="time_label", Input=nek)
tlabel.Format = "t = {time:.2f}"
tlabel_disp = Show(tlabel, rv)
tlabel_disp.WindowLocation = "Lower Right Corner"
tlabel_disp.FontSize = 28
tlabel_disp.Color = [0.1, 0.1, 0.1]

title = Text(registrationName="title")
title.Text = "PID-centered rising bubble (moving frame)"
title_disp = Show(title, rv)
title_disp.WindowLocation = "Upper Center"
title_disp.FontSize = 26
title_disp.Color = [0.1, 0.1, 0.1]

legend = Text(registrationName="legend")
legend.Text = "green: c = 0.9 in the liquid"
legend_disp = Show(legend, rv)
legend_disp.WindowLocation = "Lower Left Corner"
legend_disp.FontSize = 16
legend_disp.Color = [0.0, 0.45, 0.2]

# Fit the whole domain, seen slightly from the side and above: after the
# rotation, dolly until the domain's projected corners fill 86% of the render
# view (clear of the title and scalar bar), so taller domains are not cropped.
cx, cy = 0.5 * (bounds[0] + bounds[1]), 0.5 * (bounds[2] + bounds[3])
rv.CameraFocalPoint = [cx, cy, zmid]
rv.CameraPosition = [cx, cy, zmid + 1.0]
rv.CameraViewUp = [0.0, 1.0, 0.0]
rv.CameraViewAngle = 30.0
SetActiveView(rv)
rv.ResetCamera(False)
camera = GetActiveCamera()
camera.Azimuth(22.0)
camera.Elevation(10.0)
RENDER_FRACTION = 0.56  # render view width in the layout (SplitHorizontal below)
aspect = RENDER_FRACTION * WIDTH / HEIGHT
corners = [c + (1.0,) for c in itertools.product(bounds[0:2], bounds[2:4], bounds[4:6])]
for _ in range(10):
    M = camera.GetCompositeProjectionTransformMatrix(aspect, -1.0, 1.0)
    extent = 0.0
    for p in corners:
        x, y, z, w = M.MultiplyPoint(p)
        extent = max(extent, abs(x / w), abs(y / w))
    if abs(extent - 0.86) < 1e-3:
        break
    camera.Dolly(0.86 / extent)

# -----------------------------------------------------------------------------
# Charts: data.csv (one row per checkpoint), shown up to the current time.
# -----------------------------------------------------------------------------
csv = CSVReader(registrationName="data.csv", FileName=[CSV])

def data_up_to_time(name, columns):
    """Table of Time + columns from data.csv, with only the rows up to the
    current animation time. Only the plotted columns are passed through: the
    chart makes any column it has not seen before visible when the filter
    re-executes, so extra columns would all show up in the chart."""
    table = ProgrammableFilter(registrationName=name, Input=csv)
    table.OutputDataSetType = "vtkTable"
    # Advertise a time range so the filter re-executes whenever the animation
    # time changes.
    table.RequestInformationScript = """
from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline as sddp
oi = self.GetOutputInformation(0)
oi.Remove(sddp.TIME_STEPS())
oi.Set(sddp.TIME_RANGE(), [0.0, 1.0e6], 2)
"""
    table.Script = f"""
import numpy as np
from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline as sddp
from vtkmodules.util.numpy_support import vtk_to_numpy, numpy_to_vtk
inp = self.GetInputDataObject(0, 0)
out = self.GetOutputDataObject(0)
oi = self.GetOutputInformation(0)
t = oi.Get(sddp.UPDATE_TIME_STEP()) if oi.Has(sddp.UPDATE_TIME_STEP()) else np.inf
T = vtk_to_numpy(inp.GetColumnByName('Time'))
n = vtk_to_numpy(inp.GetColumnByName('tstep'))
# data.csv is appended to by every run in the directory. A new run (rerun or
# restart) starts where tstep or Time drops, and it supersedes the earlier rows
# after its start time Time - tstep*dt (dt from the run's own rows).
keep = np.ones(len(T), bool)
bounds = [0] + list(np.flatnonzero((np.diff(n) <= 0) | (np.diff(T) <= 0)) + 1) + [len(T)]
dt = None
for i, j in zip(bounds[:-1], bounds[1:]):
    dt = (T[j-1] - T[i]) / (n[j-1] - n[i]) if j - i > 1 else (dt or T[i] / n[i])
    keep[:i] &= T[:i] < T[i] - n[i] * dt + 0.5 * dt
# Time is written with 4 decimals, so allow for its rounding.
mask = keep & (T <= t + 5.1e-5)
Tm = T[mask]
# F_pid swings by ~+-0.015 with the level set reinit (TLSR every 0.25, CLSR
# every 0.025), which checkpoint rows 0.1 apart alias into a sawtooth. Plot its
# mean over the trailing 0.5 (two TLSR cycles) instead, computed exactly from
# the frame velocity U_frame = integral(F_pid dt).
FRAME = dict(F_pid_x="frame_u", F_pid_y="frame_v", F_pid_z="frame_w")
j = np.searchsorted(Tm, Tm - 0.51)
for name in {["Time"] + list(columns)!r}:
    v = vtk_to_numpy(inp.GetColumnByName(name))[mask]
    if name in FRAME:
        U = vtk_to_numpy(inp.GetColumnByName(FRAME[name]))[mask]
        v = np.where(j < np.arange(len(Tm)), (U - U[j]) / np.maximum(Tm - Tm[j], 1e-12), v)
    a = numpy_to_vtk(np.ascontiguousarray(v), deep=1)
    a.SetName(name)
    out.AddColumn(a)
"""
    # Keep the placeholder time range out of the animation time steps.
    tk = GetTimeKeeper()
    tk.SuppressedTimeSources = list(tk.SuppressedTimeSources) + [table]
    table.UpdatePipeline(end_time)
    return table


BLUE = ["0.12", "0.47", "0.71"]
RED = ["0.84", "0.15", "0.16"]


def series_range(table, names, pad=0.1):
    """(min, max) of the named columns of the table over the whole run, padded
    by pad times the span and always including 0, so the chart axes stay fixed
    through the animation and fit both the 2D and 3D cases."""
    from vtkmodules.util.numpy_support import vtk_to_numpy
    data = servermanager.Fetch(table)
    vals = [vtk_to_numpy(data.GetColumnByName(n)) for n in names]
    lo = min([0.0] + [float(v.min()) for v in vals if len(v)])
    hi = max([0.0] + [float(v.max()) for v in vals if len(v)])
    span = (hi - lo) or 1.0
    return lo - pad * span, hi + pad * span


def make_chart(title, ytitle, xarr_series):
    """Line chart of the given (array, label, color) series up to the current time."""
    table = data_up_to_time(f"{xarr_series[0][0]}_upto_t", [s[0] for s in xarr_series])
    yrange = series_range(table, [s[0] for s in xarr_series])
    cv = CreateView("XYChartView")
    cv.ChartTitle = title
    cv.ChartTitleFontSize = 20
    cv.LeftAxisTitle = ytitle
    cv.BottomAxisTitle = "t"
    cv.LeftAxisTitleFontSize = 16
    cv.BottomAxisTitleFontSize = 16
    cv.BottomAxisUseCustomRange = 1
    cv.BottomAxisRangeMinimum = 0.0
    cv.BottomAxisRangeMaximum = end_time
    cv.LeftAxisUseCustomRange = 1
    cv.LeftAxisRangeMinimum = yrange[0]
    cv.LeftAxisRangeMaximum = yrange[1]
    cv.LegendLocation = "TopRight"
    cv.HideTimeMarker = 0  # vertical line at the current time

    disp = Show(table, cv, "XYChartRepresentation")
    disp.UseIndexForXAxis = 0
    disp.XArrayName = "Time"
    disp.SeriesVisibility = [s[0] for s in xarr_series]
    # (paraview.simple shadows the builtin sum, so flatten with comprehensions.)
    disp.SeriesLabel = [v for s in xarr_series for v in (s[0], s[1])]
    disp.SeriesColor = [v for s in xarr_series for v in [s[0]] + s[2]]
    disp.SeriesLineThickness = [v for s in xarr_series for v in (s[0], "3")]
    return cv


force_chart = make_chart("PID centering force (per unit mass, mean over the last 0.5)", "F_pid",
        [("F_pid_x", "F_pid,x", BLUE), ("F_pid_y", "F_pid,y", RED)])
offset_chart = make_chart("Bubble centroid offset from the setpoint", "offset / D",
        [("bubble_dx", "x_c - x_0", BLUE), ("bubble_dy", "y_c - y_0", RED)])

# -----------------------------------------------------------------------------
# Layout: render view on the left, the two charts stacked on the right.
# -----------------------------------------------------------------------------
layout = CreateLayout("bubble-pid")
layout.AssignView(0, rv)
layout.SplitHorizontal(0, RENDER_FRACTION)
layout.SplitVertical(2, 0.5)
layout.AssignView(5, force_chart)
layout.AssignView(6, offset_chart)
layout.SetSize(WIDTH, HEIGHT)

scene = GetAnimationScene()
scene.UpdateAnimationUsingDataTimeSteps()
SetActiveView(rv)

# Fit the chart time axes to the loaded data whenever the animation runs, so
# the saved state also works when reloaded after the run continued (or pointed
# at another case with the same domain; the camera and glyph grid are fixed).
axis_cue = PythonAnimationCue()
axis_cue.Script = """
from paraview.simple import GetTimeKeeper, GetViews

def _fit_time_axes():
    times = GetTimeKeeper().TimestepValues
    if not times:
        return
    t_end = times[-1] if hasattr(times, '__len__') else times
    for view in GetViews():
        if view.GetXMLName() == 'XYChartView':
            view.BottomAxisRangeMaximum = t_end

def start_cue(self):
    _fit_time_axes()

def tick(self):
    _fit_time_axes()

def end_cue(self):
    pass
"""
scene.Cues.append(axis_cue)

if IN_BATCH:
    SaveState(STATE)
    print(f"saved state {STATE}")

if MAKE_MOVIE:
    if os.path.isdir(FRAMES):
        shutil.rmtree(FRAMES)
    os.makedirs(FRAMES)
    SaveAnimation(os.path.join(FRAMES, "frame.png"), layout, SaveAllViews=1,
            ImageResolution=[WIDTH, HEIGHT], FrameRate=FPS)
    print(f"saved {len(times)} frames to {FRAMES}")
    subprocess.run(["ffmpeg", "-y", "-loglevel", "error", "-framerate", str(FPS),
            "-i", os.path.join(FRAMES, "frame.%04d.png"), "-c:v", "libx264", "-crf", "18",
            "-pix_fmt", "yuv420p", MOVIE], check=True)
    print(f"saved movie {MOVIE}")
