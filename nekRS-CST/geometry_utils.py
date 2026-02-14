from dataclasses import dataclass

import gmsh

@dataclass
class XYPlaneCircle:
    center: int
    pt_12oclock: int
    pt_3oclock: int
    pt_6oclock: int
    pt_9oclock: int
    arc_12to3: int
    arc_3to6: int
    arc_6to9: int
    arc_9to12: int
    curve_loop: int
    surface: int

@dataclass
class ZCylinder:
    top: int
    bottom: int
    side: list[int]
    volume: int

@dataclass
class ZWasher:
    top: list[int]
    bottom: list[int]
    side: list[int]
    volumes: list[int]

@dataclass
class ZCylinderWithSphereHole:
    top: int
    bottom: int
    side: list[int]
    sphere: list[int]
    volumes: list[int]

g = gmsh.model.geo
m = gmsh.model.geo.mesh

def draw_xy_plane_circle(x: float, y: float, z: float, r: float, circumferential_nodes_per_quarter: int) -> XYPlaneCircle:
    center = g.add_point(x, y, z)
    pt_12oclock = g.add_point(x + r, y, z)
    pt_3oclock = g.add_point(x, y + r, z)
    pt_6oclock = g.add_point(x - r, y, z)
    pt_9oclock = g.add_point(x, y - r, z)
    arc_12to3 = g.add_circle_arc(pt_12oclock, center, pt_3oclock)
    arc_3to6 = g.add_circle_arc(pt_3oclock, center, pt_6oclock)
    arc_6to9 = g.add_circle_arc(pt_6oclock, center, pt_9oclock)
    arc_9to12 = g.add_circle_arc(pt_9oclock, center, pt_12oclock)
    g.synchronize()
    m.set_transfinite_curve(arc_12to3, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_3to6, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_6to9, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_9to12, circumferential_nodes_per_quarter)
    curve_loop = g.add_curve_loop([arc_12to3, arc_3to6, arc_6to9, arc_9to12])
    surface = g.add_plane_surface([curve_loop])
    g.synchronize()
    m.set_transfinite_surface(surface)
    m.set_recombine(2, surface)
    return XYPlaneCircle(center, pt_12oclock, pt_3oclock, pt_6oclock, pt_9oclock,
                         arc_12to3, arc_3to6, arc_6to9, arc_9to12, curve_loop, surface)

def draw_xy_plane_square(x: float, y: float, z: float, r: float, circumferential_nodes_per_quarter: int) -> XYPlaneCircle:
    center = g.add_point(x, y, z)
    pt_12oclock = g.add_point(x + r, y, z)
    pt_3oclock = g.add_point(x, y + r, z)
    pt_6oclock = g.add_point(x - r, y, z)
    pt_9oclock = g.add_point(x, y - r, z)
    arc_12to3 = g.add_line(pt_12oclock, pt_3oclock)
    arc_3to6 = g.add_line(pt_3oclock, pt_6oclock)
    arc_6to9 = g.add_line(pt_6oclock, pt_9oclock)
    arc_9to12 = g.add_line(pt_9oclock, pt_12oclock)
    g.synchronize()
    m.set_transfinite_curve(arc_12to3, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_3to6, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_6to9, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_9to12, circumferential_nodes_per_quarter)
    curve_loop = g.add_curve_loop([arc_12to3, arc_3to6, arc_6to9, arc_9to12])
    surface = g.add_plane_surface([curve_loop])
    g.synchronize()
    m.set_transfinite_surface(surface)
    m.set_recombine(2, surface)
    return XYPlaneCircle(center, pt_12oclock, pt_3oclock, pt_6oclock, pt_9oclock,
                         arc_12to3, arc_3to6, arc_6to9, arc_9to12, curve_loop, surface)

def draw_z_cylinder(top: XYPlaneCircle, bottom: XYPlaneCircle, vertical_nodes: int) -> ZCylinder:
    line_12oclock = g.add_line(top.pt_12oclock, bottom.pt_12oclock)
    line_3oclock = g.add_line(top.pt_3oclock, bottom.pt_3oclock)
    line_6oclock = g.add_line(top.pt_6oclock, bottom.pt_6oclock)
    line_9oclock = g.add_line(top.pt_9oclock, bottom.pt_9oclock)
    g.synchronize()
    m.set_transfinite_curve(line_12oclock, vertical_nodes)
    m.set_transfinite_curve(line_3oclock, vertical_nodes)
    m.set_transfinite_curve(line_6oclock, vertical_nodes)
    m.set_transfinite_curve(line_9oclock, vertical_nodes)
    curve_loop_12to3 = g.add_curve_loop([-top.arc_12to3, line_12oclock, bottom.arc_12to3, -line_3oclock])
    curve_loop_3to6 = g.add_curve_loop([-top.arc_3to6, line_3oclock, bottom.arc_3to6, -line_6oclock])
    curve_loop_6to9 = g.add_curve_loop([-top.arc_6to9, line_6oclock, bottom.arc_6to9, -line_9oclock])
    curve_loop_9to12 = g.add_curve_loop([-top.arc_9to12, line_9oclock, bottom.arc_9to12, -line_12oclock])
    surface_12to3 = g.add_surface_filling([curve_loop_12to3])
    surface_3to6 = g.add_surface_filling([curve_loop_3to6])
    surface_6to9 = g.add_surface_filling([curve_loop_6to9])
    surface_9to12 = g.add_surface_filling([curve_loop_9to12])
    g.synchronize()
    m.set_transfinite_surface(surface_12to3)
    m.set_transfinite_surface(surface_3to6)
    m.set_transfinite_surface(surface_6to9)
    m.set_transfinite_surface(surface_9to12)
    m.set_recombine(2, surface_12to3)
    m.set_recombine(2, surface_3to6)
    m.set_recombine(2, surface_6to9)
    m.set_recombine(2, surface_9to12)
    surface_loop = g.add_surface_loop([surface_12to3, surface_3to6, surface_6to9, surface_9to12, top.surface, bottom.surface])
    volume = g.add_volume([surface_loop])
    g.synchronize()
    m.set_transfinite_volume(volume)
    return ZCylinder(top.surface, bottom.surface, [surface_12to3, surface_3to6, surface_6to9, surface_9to12], volume)

def draw_z_washer(top: XYPlaneCircle, top_hole: XYPlaneCircle, bottom: XYPlaneCircle, bottom_hole: XYPlaneCircle, radial_nodes: int, vertical_nodes: int) -> ZWasher:
    line_12oclock = g.add_line(top.pt_12oclock, bottom.pt_12oclock)
    line_3oclock = g.add_line(top.pt_3oclock, bottom.pt_3oclock)
    line_6oclock = g.add_line(top.pt_6oclock, bottom.pt_6oclock)
    line_9oclock = g.add_line(top.pt_9oclock, bottom.pt_9oclock)
    g.synchronize()
    m.set_transfinite_curve(line_12oclock, vertical_nodes)
    m.set_transfinite_curve(line_3oclock, vertical_nodes)
    m.set_transfinite_curve(line_6oclock, vertical_nodes)
    m.set_transfinite_curve(line_9oclock, vertical_nodes)
    curve_loop_12to3 = g.add_curve_loop([-top.arc_12to3, line_12oclock, bottom.arc_12to3, -line_3oclock])
    curve_loop_3to6 = g.add_curve_loop([-top.arc_3to6, line_3oclock, bottom.arc_3to6, -line_6oclock])
    curve_loop_6to9 = g.add_curve_loop([-top.arc_6to9, line_6oclock, bottom.arc_6to9, -line_9oclock])
    curve_loop_9to12 = g.add_curve_loop([-top.arc_9to12, line_9oclock, bottom.arc_9to12, -line_12oclock])
    g.synchronize()
    surface_12to3 = g.add_surface_filling([curve_loop_12to3])
    surface_3to6 = g.add_surface_filling([curve_loop_3to6])
    surface_6to9 = g.add_surface_filling([curve_loop_6to9])
    surface_9to12 = g.add_surface_filling([curve_loop_9to12])
    g.synchronize()
    m.set_transfinite_surface(surface_12to3)
    m.set_transfinite_surface(surface_3to6)
    m.set_transfinite_surface(surface_6to9)
    m.set_transfinite_surface(surface_9to12)
    m.set_recombine(2, surface_12to3)
    m.set_recombine(2, surface_3to6)
    m.set_recombine(2, surface_6to9)
    m.set_recombine(2, surface_9to12)

    hole_line_12oclock = g.add_line(top_hole.pt_12oclock, bottom_hole.pt_12oclock)
    hole_line_3oclock = g.add_line(top_hole.pt_3oclock, bottom_hole.pt_3oclock)
    hole_line_6oclock = g.add_line(top_hole.pt_6oclock, bottom_hole.pt_6oclock)
    hole_line_9oclock = g.add_line(top_hole.pt_9oclock, bottom_hole.pt_9oclock)
    g.synchronize()
    m.set_transfinite_curve(hole_line_12oclock, vertical_nodes)
    m.set_transfinite_curve(hole_line_3oclock, vertical_nodes)
    m.set_transfinite_curve(hole_line_6oclock, vertical_nodes)
    m.set_transfinite_curve(hole_line_9oclock, vertical_nodes)
    hole_curve_loop_12to3 = g.add_curve_loop([-top_hole.arc_12to3, hole_line_12oclock, bottom_hole.arc_12to3, -hole_line_3oclock])
    hole_curve_loop_3to6 = g.add_curve_loop([-top_hole.arc_3to6, hole_line_3oclock, bottom_hole.arc_3to6, -hole_line_6oclock])
    hole_curve_loop_6to9 = g.add_curve_loop([-top_hole.arc_6to9, hole_line_6oclock, bottom_hole.arc_6to9, -hole_line_9oclock])
    hole_curve_loop_9to12 = g.add_curve_loop([-top_hole.arc_9to12, hole_line_9oclock, bottom_hole.arc_9to12, -hole_line_12oclock])
    hole_surface_12to3 = g.add_surface_filling([hole_curve_loop_12to3])
    hole_surface_3to6 = g.add_surface_filling([hole_curve_loop_3to6])
    hole_surface_6to9 = g.add_surface_filling([hole_curve_loop_6to9])
    hole_surface_9to12 = g.add_surface_filling([hole_curve_loop_9to12])
    g.synchronize()
    m.set_transfinite_surface(hole_surface_12to3)
    m.set_transfinite_surface(hole_surface_3to6)
    m.set_transfinite_surface(hole_surface_6to9)
    m.set_transfinite_surface(hole_surface_9to12)
    m.set_recombine(2, hole_surface_12to3)
    m.set_recombine(2, hole_surface_3to6)
    m.set_recombine(2, hole_surface_6to9)
    m.set_recombine(2, hole_surface_9to12)
    
    top_radial_12oclock = g.add_line(top_hole.pt_12oclock, top.pt_12oclock)
    top_radial_3oclock = g.add_line(top_hole.pt_3oclock, top.pt_3oclock)
    top_radial_6oclock = g.add_line(top_hole.pt_6oclock, top.pt_6oclock)
    top_radial_9oclock = g.add_line(top_hole.pt_9oclock, top.pt_9oclock)
    g.synchronize()
    m.set_transfinite_curve(top_radial_12oclock, radial_nodes)
    m.set_transfinite_curve(top_radial_3oclock, radial_nodes)
    m.set_transfinite_curve(top_radial_6oclock, radial_nodes)
    m.set_transfinite_curve(top_radial_9oclock, radial_nodes)
    top_radial_curve_loop_12to3 = g.add_curve_loop([-top_hole.arc_12to3, top_radial_12oclock, top.arc_12to3, -top_radial_3oclock])
    top_radial_curve_loop_3to6 = g.add_curve_loop([-top_hole.arc_3to6, top_radial_3oclock, top.arc_3to6, -top_radial_6oclock])
    top_radial_curve_loop_6to9 = g.add_curve_loop([-top_hole.arc_6to9, top_radial_6oclock, top.arc_6to9, -top_radial_9oclock])
    top_radial_curve_loop_9to12 = g.add_curve_loop([-top_hole.arc_9to12, top_radial_9oclock, top.arc_9to12, -top_radial_12oclock])
    top_radial_surface_12to3 = g.add_plane_surface([top_radial_curve_loop_12to3])
    top_radial_surface_3to6 = g.add_plane_surface([top_radial_curve_loop_3to6])
    top_radial_surface_6to9 = g.add_plane_surface([top_radial_curve_loop_6to9])
    top_radial_surface_9to12 = g.add_plane_surface([top_radial_curve_loop_9to12])
    g.synchronize()
    m.set_transfinite_surface(top_radial_surface_12to3)
    m.set_transfinite_surface(top_radial_surface_3to6)
    m.set_transfinite_surface(top_radial_surface_6to9)
    m.set_transfinite_surface(top_radial_surface_9to12)
    m.set_recombine(2, top_radial_surface_12to3)
    m.set_recombine(2, top_radial_surface_3to6)
    m.set_recombine(2, top_radial_surface_6to9)
    m.set_recombine(2, top_radial_surface_9to12)

    bottom_radial_12oclock = g.add_line(bottom_hole.pt_12oclock, bottom.pt_12oclock)
    bottom_radial_3oclock = g.add_line(bottom_hole.pt_3oclock, bottom.pt_3oclock)
    bottom_radial_6oclock = g.add_line(bottom_hole.pt_6oclock, bottom.pt_6oclock)
    bottom_radial_9oclock = g.add_line(bottom_hole.pt_9oclock, bottom.pt_9oclock)
    g.synchronize()
    m.set_transfinite_curve(bottom_radial_12oclock, radial_nodes)
    m.set_transfinite_curve(bottom_radial_3oclock, radial_nodes)
    m.set_transfinite_curve(bottom_radial_6oclock, radial_nodes)
    m.set_transfinite_curve(bottom_radial_9oclock, radial_nodes)
    bottom_radial_curve_loop_12to3 = g.add_curve_loop([-bottom_hole.arc_12to3, bottom_radial_12oclock, bottom.arc_12to3, -bottom_radial_3oclock])
    bottom_radial_curve_loop_3to6 = g.add_curve_loop([-bottom_hole.arc_3to6, bottom_radial_3oclock, bottom.arc_3to6, -bottom_radial_6oclock])
    bottom_radial_curve_loop_6to9 = g.add_curve_loop([-bottom_hole.arc_6to9, bottom_radial_6oclock, bottom.arc_6to9, -bottom_radial_9oclock])
    bottom_radial_curve_loop_9to12 = g.add_curve_loop([-bottom_hole.arc_9to12, bottom_radial_9oclock, bottom.arc_9to12, -bottom_radial_12oclock])
    bottom_radial_surface_12to3 = g.add_plane_surface([bottom_radial_curve_loop_12to3])
    bottom_radial_surface_3to6 = g.add_plane_surface([bottom_radial_curve_loop_3to6])
    bottom_radial_surface_6to9 = g.add_plane_surface([bottom_radial_curve_loop_6to9])
    bottom_radial_surface_9to12 = g.add_plane_surface([bottom_radial_curve_loop_9to12])
    g.synchronize()
    m.set_transfinite_surface(bottom_radial_surface_12to3)
    m.set_transfinite_surface(bottom_radial_surface_3to6)
    m.set_transfinite_surface(bottom_radial_surface_6to9)
    m.set_transfinite_surface(bottom_radial_surface_9to12)
    m.set_recombine(2, bottom_radial_surface_12to3)
    m.set_recombine(2, bottom_radial_surface_3to6)
    m.set_recombine(2, bottom_radial_surface_6to9)
    m.set_recombine(2, bottom_radial_surface_9to12)

    separator_curve_loop_12 = g.add_curve_loop([-line_12oclock, -top_radial_12oclock, hole_line_12oclock, bottom_radial_12oclock])
    separator_curve_loop_3 = g.add_curve_loop([-line_3oclock, -top_radial_3oclock, hole_line_3oclock, bottom_radial_3oclock])
    separator_curve_loop_6 = g.add_curve_loop([-line_6oclock, -top_radial_6oclock, hole_line_6oclock, bottom_radial_6oclock])
    separator_curve_loop_9 = g.add_curve_loop([-line_9oclock, -top_radial_9oclock, hole_line_9oclock, bottom_radial_9oclock])
    separator_12 = g.add_plane_surface([separator_curve_loop_12])
    separator_3 = g.add_plane_surface([separator_curve_loop_3])
    separator_6 = g.add_plane_surface([separator_curve_loop_6])
    separator_9 = g.add_plane_surface([separator_curve_loop_9])
    g.synchronize()
    m.set_transfinite_surface(separator_12)
    m.set_transfinite_surface(separator_3)
    m.set_transfinite_surface(separator_6)
    m.set_transfinite_surface(separator_9)
    m.set_recombine(2, separator_12)
    m.set_recombine(2, separator_3)
    m.set_recombine(2, separator_6)
    m.set_recombine(2, separator_9)

    surface_loop_12to3 = g.add_surface_loop([surface_12to3, hole_surface_12to3, top_radial_surface_12to3, bottom_radial_surface_12to3, separator_12, separator_3])
    surface_loop_3to6 = g.add_surface_loop([surface_3to6, hole_surface_3to6, top_radial_surface_3to6, bottom_radial_surface_3to6, separator_3, separator_6])
    surface_loop_6to9 = g.add_surface_loop([surface_6to9, hole_surface_6to9, top_radial_surface_6to9, bottom_radial_surface_6to9, separator_6, separator_9])
    surface_loop_9to12 = g.add_surface_loop([surface_9to12, hole_surface_9to12, top_radial_surface_9to12, bottom_radial_surface_9to12, separator_9, separator_12])
    volume_12to3 = g.add_volume([surface_loop_12to3])
    volume_3to6 = g.add_volume([surface_loop_3to6])
    volume_6to9 = g.add_volume([surface_loop_6to9])
    volume_9to12 = g.add_volume([surface_loop_9to12])
    g.synchronize()
    m.set_transfinite_volume(volume_12to3)
    m.set_transfinite_volume(volume_3to6)
    m.set_transfinite_volume(volume_6to9)
    m.set_transfinite_volume(volume_9to12)

    outer_surfaces = [surface_12to3, surface_3to6, surface_6to9, surface_9to12]
    top_radial_surfaces = [top_radial_surface_12to3, top_radial_surface_3to6, top_radial_surface_6to9, top_radial_surface_9to12]
    bottom_radial_surfaces = [bottom_radial_surface_12to3, bottom_radial_surface_3to6, bottom_radial_surface_6to9, bottom_radial_surface_9to12]
    volumes = [volume_12to3, volume_3to6, volume_6to9, volume_9to12]

    return ZWasher(top_radial_surfaces, bottom_radial_surfaces, outer_surfaces, volumes)

def draw_z_cylinder_with_sphere_hole(top: XYPlaneCircle, bottom: XYPlaneCircle, sphere_radius: float, circumferential_nodes_per_quarter: int, vertical_nodes: int) -> ZCylinderWithSphereHole:
    """
    Create a cylinder with a spherical hole at the center.
    Returns a ZCylinderWithSphereHole dataclass.
    """
    # Cylinder sides
    cyl_line_12oclock = g.add_line(top.pt_12oclock, bottom.pt_12oclock)
    cyl_line_3oclock = g.add_line(top.pt_3oclock, bottom.pt_3oclock)
    cyl_line_6oclock = g.add_line(top.pt_6oclock, bottom.pt_6oclock)
    cyl_line_9oclock = g.add_line(top.pt_9oclock, bottom.pt_9oclock)
    g.synchronize()
    m.set_transfinite_curve(cyl_line_12oclock, vertical_nodes)
    m.set_transfinite_curve(cyl_line_3oclock, vertical_nodes)
    m.set_transfinite_curve(cyl_line_6oclock, vertical_nodes)
    m.set_transfinite_curve(cyl_line_9oclock, vertical_nodes)
    cyl_side_12to3 = g.add_surface_filling([g.add_curve_loop([-top.arc_12to3, cyl_line_12oclock, bottom.arc_12to3, -cyl_line_3oclock])])
    cyl_side_3to6 = g.add_surface_filling([g.add_curve_loop([-top.arc_3to6, cyl_line_3oclock, bottom.arc_3to6, -cyl_line_6oclock])])
    cyl_side_6to9 = g.add_surface_filling([g.add_curve_loop([-top.arc_6to9, cyl_line_6oclock, bottom.arc_6to9, -cyl_line_9oclock])])
    cyl_side_9to12 = g.add_surface_filling([g.add_curve_loop([-top.arc_9to12, cyl_line_9oclock, bottom.arc_9to12, -cyl_line_12oclock])])
    g.synchronize()
    m.set_transfinite_surface(cyl_side_12to3)
    m.set_transfinite_surface(cyl_side_3to6)
    m.set_transfinite_surface(cyl_side_6to9)
    m.set_transfinite_surface(cyl_side_9to12)
    m.set_recombine(2, cyl_side_12to3)
    m.set_recombine(2, cyl_side_3to6)
    m.set_recombine(2, cyl_side_6to9)
    m.set_recombine(2, cyl_side_9to12)

    # Sphere points (centered at origin)
    x = sphere_radius * (2 ** 0.5) / 2
    center = g.add_point(0, 0, 0)
    top_12 = g.add_point(x, 0, x)
    top_3 = g.add_point(0, x, x)
    top_6 = g.add_point(-x, 0, x)
    top_9 = g.add_point(0, -x, x)
    bottom_12 = g.add_point(x, 0, -x)
    bottom_3 = g.add_point(0, x, -x)
    bottom_6 = g.add_point(-x, 0, -x)
    bottom_9 = g.add_point(0, -x, -x)

    # Sphere top and bottom
    arc_top_12to3 = g.add_circle_arc(top_12, center, top_3)
    arc_top_3to6 = g.add_circle_arc(top_3, center, top_6)
    arc_top_6to9 = g.add_circle_arc(top_6, center, top_9)
    arc_top_9to12 = g.add_circle_arc(top_9, center, top_12)
    g.synchronize()
    m.set_transfinite_curve(arc_top_12to3, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_top_3to6, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_top_6to9, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_top_9to12, circumferential_nodes_per_quarter)
    sphere_top = g.add_surface_filling([g.add_curve_loop([arc_top_12to3, arc_top_3to6, arc_top_6to9, arc_top_9to12])], sphereCenterTag=center)
    g.synchronize()
    m.set_transfinite_surface(sphere_top)
    m.set_recombine(2, sphere_top)
    arc_bot_12to3 = g.add_circle_arc(bottom_12, center, bottom_3)
    arc_bot_3to6 = g.add_circle_arc(bottom_3, center, bottom_6)
    arc_bot_6to9 = g.add_circle_arc(bottom_6, center, bottom_9)
    arc_bot_9to12 = g.add_circle_arc(bottom_9, center, bottom_12)
    g.synchronize()
    m.set_transfinite_curve(arc_bot_12to3, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_bot_3to6, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_bot_6to9, circumferential_nodes_per_quarter)
    m.set_transfinite_curve(arc_bot_9to12, circumferential_nodes_per_quarter)
    sphere_bottom = g.add_surface_filling([g.add_curve_loop([arc_bot_12to3, arc_bot_3to6, arc_bot_6to9, arc_bot_9to12])], sphereCenterTag=center)
    g.synchronize()
    m.set_transfinite_surface(sphere_bottom)
    m.set_recombine(2, sphere_bottom)

    # Sphere sides
    arc_top_3_bot_3 = g.add_circle_arc(top_3, center, bottom_3)
    arc_top_12_bot_12 = g.add_circle_arc(top_12, center, bottom_12)
    arc_top_6_bot_6 = g.add_circle_arc(top_6, center, bottom_6)
    arc_top_9_bot_9 = g.add_circle_arc(top_9, center, bottom_9)
    g.synchronize()
    m.set_transfinite_curve(arc_top_3_bot_3, vertical_nodes)
    m.set_transfinite_curve(arc_top_12_bot_12, vertical_nodes)
    m.set_transfinite_curve(arc_top_6_bot_6, vertical_nodes)
    m.set_transfinite_curve(arc_top_9_bot_9, vertical_nodes)
    sphere_12to3 = g.add_surface_filling([g.add_curve_loop([arc_top_12to3, arc_top_3_bot_3, -arc_bot_12to3, -arc_top_12_bot_12])], sphereCenterTag=center)
    g.synchronize()
    m.set_transfinite_surface(sphere_12to3)
    m.set_recombine(2, sphere_12to3)
    sphere_3to6 = g.add_surface_filling([g.add_curve_loop([arc_top_3to6, arc_top_6_bot_6, -arc_bot_3to6, -arc_top_3_bot_3])], sphereCenterTag=center)
    g.synchronize()
    m.set_transfinite_surface(sphere_3to6)
    m.set_recombine(2, sphere_3to6)
    sphere_6to9 = g.add_surface_filling([g.add_curve_loop([arc_top_6to9, arc_top_9_bot_9, -arc_bot_6to9, -arc_top_6_bot_6])], sphereCenterTag=center)
    g.synchronize()
    m.set_transfinite_surface(sphere_6to9)
    m.set_recombine(2, sphere_6to9)
    sphere_9to12 = g.add_surface_filling([g.add_curve_loop([arc_top_9to12, arc_top_12_bot_12, -arc_bot_9to12, -arc_top_9_bot_9])], sphereCenterTag=center)
    g.synchronize()
    m.set_transfinite_surface(sphere_9to12)
    m.set_recombine(2, sphere_9to12)

    # Separator surfaces
    sphere_to_cylinder_12_top = g.add_line(top_12, top.pt_12oclock)
    sphere_to_cylinder_3_top = g.add_line(top_3, top.pt_3oclock)
    sphere_to_cylinder_6_top = g.add_line(top_6, top.pt_6oclock)
    sphere_to_cylinder_9_top = g.add_line(top_9, top.pt_9oclock)
    sphere_to_cylinder_12_bot = g.add_line(bottom_12, bottom.pt_12oclock)
    sphere_to_cylinder_3_bot = g.add_line(bottom_3, bottom.pt_3oclock)
    sphere_to_cylinder_6_bot = g.add_line(bottom_6, bottom.pt_6oclock)
    sphere_to_cylinder_9_bot = g.add_line(bottom_9, bottom.pt_9oclock)
    g.synchronize()
    m.set_transfinite_curve(sphere_to_cylinder_12_top, vertical_nodes)
    m.set_transfinite_curve(sphere_to_cylinder_3_top, vertical_nodes)
    m.set_transfinite_curve(sphere_to_cylinder_6_top, vertical_nodes)
    m.set_transfinite_curve(sphere_to_cylinder_9_top, vertical_nodes)
    m.set_transfinite_curve(sphere_to_cylinder_12_bot, vertical_nodes)
    m.set_transfinite_curve(sphere_to_cylinder_3_bot, vertical_nodes)
    m.set_transfinite_curve(sphere_to_cylinder_6_bot, vertical_nodes)
    m.set_transfinite_curve(sphere_to_cylinder_9_bot, vertical_nodes)
    separator_top_12to3 = g.add_plane_surface([g.add_curve_loop([-sphere_to_cylinder_12_top, arc_top_12to3, sphere_to_cylinder_3_top, -top.arc_12to3])])
    separator_top_3to6 = g.add_plane_surface([g.add_curve_loop([-sphere_to_cylinder_3_top, arc_top_3to6, sphere_to_cylinder_6_top, -top.arc_3to6])])
    separator_top_6to9 = g.add_plane_surface([g.add_curve_loop([-sphere_to_cylinder_6_top, arc_top_6to9, sphere_to_cylinder_9_top, -top.arc_6to9])])
    separator_top_9to12 = g.add_plane_surface([g.add_curve_loop([-sphere_to_cylinder_9_top, arc_top_9to12, sphere_to_cylinder_12_top, -top.arc_9to12])])
    separator_bot_12to3 = g.add_plane_surface([g.add_curve_loop([-sphere_to_cylinder_12_bot, arc_bot_12to3, sphere_to_cylinder_3_bot, -bottom.arc_12to3])])
    separator_bot_3to6 = g.add_plane_surface([g.add_curve_loop([-sphere_to_cylinder_3_bot, arc_bot_3to6, sphere_to_cylinder_6_bot, -bottom.arc_3to6])])
    separator_bot_6to9 = g.add_plane_surface([g.add_curve_loop([-sphere_to_cylinder_6_bot, arc_bot_6to9, sphere_to_cylinder_9_bot, -bottom.arc_6to9])])
    separator_bot_9to12 = g.add_plane_surface([g.add_curve_loop([-sphere_to_cylinder_9_bot, arc_bot_9to12, sphere_to_cylinder_12_bot, -bottom.arc_9to12])])
    separator_12 = g.add_plane_surface([g.add_curve_loop([-cyl_line_12oclock, -sphere_to_cylinder_12_top, arc_top_12_bot_12, sphere_to_cylinder_12_bot])])
    separator_3 = g.add_plane_surface([g.add_curve_loop([-cyl_line_3oclock, -sphere_to_cylinder_3_top, arc_top_3_bot_3, sphere_to_cylinder_3_bot])])
    separator_6 = g.add_plane_surface([g.add_curve_loop([-cyl_line_6oclock, -sphere_to_cylinder_6_top, arc_top_6_bot_6, sphere_to_cylinder_6_bot])])
    separator_9 = g.add_plane_surface([g.add_curve_loop([-cyl_line_9oclock, -sphere_to_cylinder_9_top, arc_top_9_bot_9, sphere_to_cylinder_9_bot])])
    g.synchronize()
    separators = [separator_top_12to3, separator_top_3to6, separator_top_6to9, separator_top_9to12,
                  separator_bot_12to3, separator_bot_3to6, separator_bot_6to9, separator_bot_9to12,
                  separator_12, separator_3, separator_6, separator_9]
    for s in separators:
        m.set_transfinite_surface(s)
        m.set_recombine(2, s)

    # Volumes
    volume_top = g.add_volume([g.add_surface_loop([sphere_top, separator_top_12to3, separator_top_3to6, separator_top_6to9, separator_top_9to12, top.surface])])
    volume_bottom = g.add_volume([g.add_surface_loop([sphere_bottom, separator_bot_12to3, separator_bot_3to6, separator_bot_6to9, separator_bot_9to12, bottom.surface])])
    volume_12to3 = g.add_volume([g.add_surface_loop([cyl_side_12to3, sphere_12to3, separator_top_12to3, separator_bot_12to3, separator_12, separator_3])])
    volume_3to6 = g.add_volume([g.add_surface_loop([cyl_side_3to6, sphere_3to6, separator_top_3to6, separator_bot_3to6, separator_3, separator_6])])
    volume_6to9 = g.add_volume([g.add_surface_loop([cyl_side_6to9, sphere_6to9, separator_top_6to9, separator_bot_6to9, separator_6, separator_9])])
    volume_9to12 = g.add_volume([g.add_surface_loop([cyl_side_9to12, sphere_9to12, separator_top_9to12, separator_bot_9to12, separator_9, separator_12])])
    g.synchronize()
    volumes = [volume_top, volume_bottom, volume_12to3, volume_3to6, volume_6to9, volume_9to12]
    for v in volumes:
        m.set_transfinite_volume(v)

    return ZCylinderWithSphereHole(top.surface, bottom.surface,
                                   [cyl_side_12to3, cyl_side_3to6, cyl_side_6to9, cyl_side_9to12],
                                   [sphere_top, sphere_bottom, sphere_12to3, sphere_3to6, sphere_6to9, sphere_9to12],
                                   volumes)