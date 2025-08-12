// gmsh geometry file for a cylindrical domain with a spherical cavity,
// representing bulk liquid with a bubble void centered at (0, 0).
//
// We take the z-axis to be the vertical direction, and flow to be
// bottom-to-top as our convention (to match the original .exo mesh). Units
// are currently arbitrary; this file will be updated to use SI units when we
// input the parameters for our desired physical system.
//
// For now we have uniform meshing throughout our domain; this file will be
// later updated to refine the region adjacent to the bubble.

// Parameters.
bubble_radius = 0.5;
domain_radius = 10;
domain_upstream_height = 10;
domain_downstream_height = 20;

// Circle for the liquid inlet.
Point(11) = {0, 0, -domain_upstream_height}; // center
Point(12) = {domain_radius, 0, -domain_upstream_height}; // 12 o'clock
Point(13) = {0, domain_radius, -domain_upstream_height}; // 3 o'clock
Point(14) = {-domain_radius, 0, -domain_upstream_height}; // 6 o'clock
Point(15) = {0, -domain_radius, -domain_upstream_height}; // 9 o'clock
// The Circle primitive doesn't support constructing a full circle in one go...
// need to build it in quarters.
Circle(101) = {12, 11, 13};
Circle(102) = {13, 11, 14};
Circle(103) = {14, 11, 15};
Circle(104) = {15, 11, 12};
Curve Loop(105) = {101, 102, 103, 104};
// Make the final Surface.
Surface(1001) = {105};

// Circle for the liquid outlet.
Point(21) = {0, 0, domain_downstream_height}; // center
Point(22) = {domain_radius, 0, domain_downstream_height}; // 12 o'clock
Point(23) = {0, domain_radius, domain_downstream_height}; // 3 o'clock
Point(24) = {-domain_radius, 0, domain_downstream_height}; // 6 o'clock
Point(25) = {0, -domain_radius, domain_downstream_height}; // 9 o'clock
// The Circle primitive doesn't support constructing a full circle in one go...
// need to build it in quarters.
Circle(201) = {22, 21, 23};
Circle(202) = {23, 21, 24};
Circle(203) = {24, 21, 25};
Circle(204) = {25, 21, 22};
Curve Loop(205) = {201, 202, 203, 204};
// Make the final Surface.
Surface(1002) = {205};

// Side wall of the cylindrical domain.
Surface(1003) = {105, 205};

// Bubble surface.
Point(31) = {0, 0, 0}; // center
Point(32) = {bubble_radius, 0, 0}; // 12 o'clock
Point(33) = {0, bubble_radius, 0, 0}; // 3 o'clock
Point(34) = {-bubble_radius, 0, 0}; // 6 o'clock
Point(35) = {0, -bubble_radius, 0, 0}; // 9 o'clock
// The Circle primitive doesn't support constructing a full circle in one go...
// need to build it in quarters.
Circle(301) = {32, 31, 33};
Circle(302) = {33, 31, 34};
Circle(303) = {34, 31, 35};
Circle(304) = {35, 31, 32};
Curve Loop(305) = {301, 302, 303, 304};
Surface(1004) = {305} In Sphere {31};

// Volume to mesh on.
Surface Loop(1005) = {1001, 1003, 1002};
Volume(10001) = {1005, 1004};

// Physical surfaces and volumes, for meshing and export to OpenFOAM.
Physical Surface("inlet") = {1001};
Physical Surface("outlet") = {1002};
Physical Surface("wall") = {1003};
Physical Surface("bubble") = {1004};
Physical Volume("fluid") = {10001};

// Set mesh granularity.
Mesh.CharacteristicLengthFactor = 0.1;
// In the future, characteristic length can be set in different volumes for
// focused refinement.
//Characteristic Length{ PointsOf{ Volume{10001}; } } = 0.1;
