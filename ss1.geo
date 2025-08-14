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
//
//  $ gmsh -3 -order 2 -format msh2
//  $ gmsh2nek

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

// Side wall. We make it in 4 parts so that the axis dimension is properly
// subdivided in 1D.
Line(301) = {12, 22};
Line(302) = {13, 23};
Line(303) = {14, 24};
Line(304) = {15, 25};
Curve Loop(305) = {301, 201, -302, -101};
Curve Loop(306) = {302, 202, -303, -102};
Curve Loop(307) = {303, 203, -304, -103};
Curve Loop(308) = {304, 204, -301, -104};
Surface(1003) = {305};
Surface(1004) = {306};
Surface(1005) = {307};
Surface(1006) = {308};

// Bubble sphere surface.
// Making a sphere is kind of a chore using only builtins (no OpenCASCADE).
Point(41) = {0, 0, 0}; // center
Point(42) = {bubble_radius, 0, 0}; // 12 o'clock
Point(43) = {-bubble_radius, 0, 0}; // 6 o'clock
Point(44) = {0, bubble_radius, 0, 0}; // 3 o'clock
Point(45) = {0, -bubble_radius, 0, 0}; // 9 o'clock
Point(46) = {0, 0, bubble_radius}; // north pole
Point(47) = {0, 0, -bubble_radius}; // south pole
Circle(401) = {42, 41, 44};
Circle(402) = {44, 41, 43};
Circle(403) = {43, 41, 45};
Circle(404) = {45, 41, 42};
Circle(405) = {42, 41, 46};
Circle(406) = {46, 41, 43};
Circle(407) = {43, 41, 47};
Circle(408) = {47, 41, 42};
Circle(409) = {44, 41, 46};
Circle(410) = {46, 41, 45};
Circle(411) = {45, 41, 47};
Circle(412) = {47, 41, 44};
Curve Loop(413) = {402, 407, 412};
Surface(1007) = {413} In Sphere {41};
Curve Loop(414) = {-402, 409, 406};
Surface(1008) = {414} In Sphere {41};
Curve Loop(415) = {-403, -406, 410};
Surface(1009) = {415} In Sphere {41};
Curve Loop(416) = {-403, 407, -411};
Surface(1010) = {416} In Sphere {41};
Curve Loop(417) = {-404, 411, 408};
Surface(1011) = {417} In Sphere {41};
Curve Loop(418) = {404, 405, 410};
Surface(1012) = {418} In Sphere {41};
Curve Loop(419) = {-401, 405, -409};
Surface(1013) = {419} In Sphere {41};
Curve Loop(420) = {-401, -408, 412};
Surface(1014) = {420} In Sphere {41};

// Volume to be meshed.
Surface Loop(1015) = {1001, 1002, 1003, 1004, 1005, 1006}; // outer surface
Surface Loop(1016) = {1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014}; // hole
Volume(10001) = {1015, 1016};

// Recombine tris into quads. (This just generates a pyramid/tet hybrid mesh
// though.)
//Transfinite Surface{1001, 1002, 1003, 1004, 1005, 1006};
//Recombine Surface{1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014};

// Physical surfaces and volumes, for meshing and export to OpenFOAM.
Physical Surface("inlet") = {1001};
Physical Surface("outlet") = {1002};
Physical Surface("wall") = {1003, 1004, 1005, 1006};
Physical Surface("bubble") = {1007, 1008, 1009, 1010, 1011, 1012, 1013, 1014};
Physical Volume("fluid") = {10001};

// Define a simple line for sphere wake (used for length scaling).
Point(51) = {0, 0, 2*bubble_radius};
Curve(501) = {41, 51};

// Set length scale field based on distance from sphere center or sphere wake
// line.
// Distance of bubble_radius => length of 0.01.
// Distance of 1.2*bubble_radius => length of 0.15.
Field[1] = Distance;
Field[1].CurvesList = {501};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.01;
Field[2].SizeMax = 0.15;
Field[2].DistMin = bubble_radius;
Field[2].DistMax = bubble_radius*1.2;
Background Field = 2;

// Global refinement scale. (Divide by 2 for real scale, due to hex27/order 2
// elements.)
Mesh.CharacteristicLengthFactor = 10;

// Ignore curvature etc -- use purely our length scale field.
Mesh.MeshSizeFromPoints = 0;
Mesh.MeshSizeFromCurvature = 0;
Mesh.MeshSizeExtendFromBoundary = 0;
