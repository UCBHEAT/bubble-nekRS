import gmsh

from geometry_utils import (
    draw_xy_plane_circle, draw_xy_plane_square, draw_z_cylinder, draw_z_washer, draw_z_cylinder_with_sphere_hole
)

# Parameters
bubble_radius = 0.5
refined_zone_radius = 1.5
domain_radius = 10
domain_downstream_height = 20
domain_upstream_height = 10  # upstream is on the bottom

node_density_coarse = 2
node_density_fine = 4
node_density_ultrafine = 6

vertical_nodes_1to2 = int(node_density_coarse*(domain_downstream_height*3/4))
vertical_nodes_2to3 = int(node_density_ultrafine*(domain_downstream_height/4 - refined_zone_radius))
vertical_nodes_3to4 = int(node_density_fine*(2*refined_zone_radius))
vertical_nodes_4to5 = int(node_density_coarse*(domain_upstream_height - refined_zone_radius))
radial_nodes_outer = int(node_density_coarse*(domain_radius - refined_zone_radius))
radial_nodes_inner = int(node_density_fine*refined_zone_radius)
circumferential_nodes_per_quarter = 30

gmsh.initialize()
gmsh.model.add("ss1")

# Layer 1
layer_1_circle = draw_xy_plane_circle(0, 0, domain_downstream_height, domain_radius, circumferential_nodes_per_quarter)
layer_1_hole = draw_xy_plane_square(0, 0, domain_downstream_height, refined_zone_radius, circumferential_nodes_per_quarter)

# Layer 2
layer_2_circle = draw_xy_plane_circle(0, 0, domain_downstream_height/4, domain_radius, circumferential_nodes_per_quarter)
layer_2_hole = draw_xy_plane_square(0, 0, domain_downstream_height/4, refined_zone_radius, circumferential_nodes_per_quarter)

# Layer 3
washer_top = draw_xy_plane_circle(0, 0, refined_zone_radius, domain_radius, circumferential_nodes_per_quarter)
washer_hole_top = draw_xy_plane_square(0, 0, refined_zone_radius, refined_zone_radius, circumferential_nodes_per_quarter)

# Layer 4
washer_bottom = draw_xy_plane_circle(0, 0, -refined_zone_radius, domain_radius, circumferential_nodes_per_quarter)
washer_hole_bottom = draw_xy_plane_square(0, 0, -refined_zone_radius, refined_zone_radius, circumferential_nodes_per_quarter)

# Layer 5
layer_5_circle = draw_xy_plane_circle(0, 0, -domain_upstream_height, domain_radius, circumferential_nodes_per_quarter)
layer_5_hole = draw_xy_plane_square(0, 0, -domain_upstream_height, refined_zone_radius, circumferential_nodes_per_quarter)

# Volumes
layer_1to2_outer = draw_z_washer(layer_1_circle, layer_1_hole, layer_2_circle, layer_2_hole, radial_nodes_outer, vertical_nodes_1to2)
layer_1to2_inner = draw_z_cylinder(layer_1_hole, layer_2_hole, vertical_nodes_1to2)
layer_2to3_outer = draw_z_washer(layer_2_circle, layer_2_hole, washer_top, washer_hole_top, radial_nodes_outer, vertical_nodes_2to3)
layer_2to3_inner = draw_z_cylinder(layer_2_hole, washer_hole_top, vertical_nodes_2to3)
layer_3to4_outer = draw_z_washer(washer_top, washer_hole_top, washer_bottom, washer_hole_bottom, radial_nodes_outer, vertical_nodes_3to4)
layer_3to4_inner = draw_z_cylinder_with_sphere_hole(washer_hole_top, washer_hole_bottom, bubble_radius, circumferential_nodes_per_quarter, vertical_nodes_3to4)
layer_4to5_outer = draw_z_washer(washer_bottom, washer_hole_bottom, layer_5_circle, layer_5_hole, radial_nodes_outer, vertical_nodes_4to5)
layer_4to5_inner = draw_z_cylinder(washer_hole_bottom, layer_5_hole, vertical_nodes_4to5)

gmsh.model.geo.synchronize()
gmsh.model.addPhysicalGroup(2, [layer_1to2_inner.top] + layer_1to2_outer.top, name="inlet")
gmsh.model.addPhysicalGroup(2, [layer_4to5_inner.bottom] + layer_4to5_outer.bottom, name="outlet")
gmsh.model.addPhysicalGroup(2, layer_1to2_inner.side + layer_2to3_outer.side + layer_3to4_outer.side + layer_4to5_outer.side, name="walls")
gmsh.model.addPhysicalGroup(2, layer_3to4_inner.sphere, name="sphere")
gmsh.model.addPhysicalGroup(3, [*layer_1to2_outer.volumes, layer_1to2_inner.volume, *layer_2to3_outer.volumes,
                                layer_2to3_inner.volume, *layer_3to4_outer.volumes, *layer_3to4_inner.volumes,
                                *layer_4to5_outer.volumes, layer_4to5_inner.volume], name="fluid")

# Uncomment to inspect generated objects, for debugging meshing issues
#gmsh.write("ss1.geo_unrolled")

gmsh.model.mesh.generate(3)
gmsh.model.mesh.set_order(2)
gmsh.option.setNumber("Mesh.MshFileVersion", 2.2)
gmsh.write("ss1.msh")

# Uncomment to visualize geometry/mesh in GUI
#gmsh.fltk.run()

gmsh.finalize()
