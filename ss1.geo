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

// Side wall. We make it in 4 parts so that they can be part of transfinite
// volumes with the sphere pieces.
Line(301) = {12, 22};
Line(302) = {13, 23};
Line(303) = {14, 24};
Line(304) = {15, 25};
Curve Loop(305) = {301, 201, -302, -101};
Curve Loop(306) = {302, 202, -303, -102};
Curve Loop(307) = {303, 203, -304, -103};
Curve Loop(308) = {304, 204, -301, -104};
Surface(3001) = {305};
Surface(3002) = {306};
Surface(3003) = {307};
Surface(3004) = {308};

// Bubble sphere surface.
// Making a sphere is kind of a chore using only builtins (no OpenCASCADE).
x = bubble_radius * Sqrt(2)/2;
Point(41) = {0, 0, 0}; // center
Point(42) = {x, 0, x}; // top 12 o'clock
Point(43) = {0, x, x}; // top 3 o'clock
Point(44) = {-x, 0, x}; // top 6 o'clock
Point(45) = {0, -x, x}; // top 9 o'clock
Point(46) = {x, 0, -x}; // bottom 12 o'clock
Point(47) = {0, x, -x}; // bottom 3 o'clock
Point(48) = {-x, 0, -x}; // bottom 6 o'clock
Point(49) = {0, -x, -x}; // bottom 9 o'clock
// top face
Circle(401) = {42, 41, 43};
Circle(402) = {43, 41, 44};
Circle(403) = {44, 41, 45};
Circle(404) = {45, 41, 42};
Curve Loop(405) = {401, 402, 403, 404};
Surface(4001) = {405};
// bottom face
Circle(426) = {46, 41, 47};
Circle(427) = {47, 41, 48};
Circle(428) = {48, 41, 49};
Circle(429) = {49, 41, 46};
Curve Loop(430) = {426, 427, 428, 429};
Surface(4006) = {430};
// 12-3 o'clock
//Circle(406) = {42, 41, 43};
Circle(407) = {43, 41, 47};
//Circle(408) = {47, 41, 46};
//Circle(409) = {46, 41, 42};
Circle(422) = {42, 41, 46};
Curve Loop(410) = {401, 407, -426, -422};
Surface(4002) = {410};
// 3-6 o'clock
//Circle(411) = {43, 41, 44};
Circle(412) = {44, 41, 48};
//Circle(413) = {48, 41, 47};
//Circle(414) = {47, 41, 43};
Curve Loop(415) = {402, 412, -427, -407};
Surface(4003) = {415};
// 6-9 o'clock
//Circle(416) = {44, 41, 45};
Circle(417) = {45, 41, 49};
//Circle(418) = {49, 41, 48};
//Circle(419) = {48, 41, 44};
Curve Loop(420) = {403, 417, -428, -412};
Surface(4004) = {420};
// 9-12 o'clock
//Circle(421) = {45, 41, 42};
//Circle(422) = {42, 41, 46};
//Circle(423) = {46, 41, 49};
//Circle(424) = {49, 41, 45};
Curve Loop(425) = {404, 422, -429, -417};
Surface(4005) = {425};

// Top connector lines & surfaces
Line(501) = {22, 42};
Line(502) = {23, 43};
Line(503) = {24, 44};
Line(504) = {25, 45};
Curve Loop(509) = {201, 502, -401, -501};
Surface(5001) = {509};
Curve Loop(510) = {202, 503, -402, -502};
Surface(5002) = {510};
Curve Loop(511) = {203, 504, -403, -503};
Surface(5003) = {511};
Curve Loop(512) = {204, 501, -404, -504};
Surface(5004) = {512};

// Bottom connector lines & surfaces
Line(505) = {12, 46};
Line(506) = {13, 47};
Line(507) = {14, 48};
Line(508) = {15, 49};
Curve Loop(513) = {101, 506, -426, -505};
Surface(5005) = {513};
Curve Loop(514) = {102, 507, -427, -506};
Surface(5006) = {514};
Curve Loop(515) = {103, 508, -428, -507};
Surface(5007) = {515};
Curve Loop(516) = {104, 505, -429, -508};
Surface(5008) = {516};

// Side wedges
Curve Loop(517) = {501, 422, -505, 301}; // 12 o'clock
Surface(5009) = {517};
Curve Loop(518) = {502, 407, -506, 302}; // 3 o'clock
Surface(5010) = {518};
Curve Loop(519) = {503, 412, -507, 303}; // 6 o'clock
Surface(5011) = {519};
Curve Loop(520) = {504, 417, -508, 304}; // 9 o'clock
Surface(5012) = {520};

// Volumes
// Top
Surface Loop(6001) = {1002, 5001, 5002, 5003, 5004, 4001};
Volume(60001) = {6001};
// Bottom
Surface Loop(6002) = {5007, 4006, 5005, 5006, 5008, 1001};
Volume(60002) = {6002};
// 12-3 o'clock
Surface Loop(6003) = {5009, 5010, 4002, 3001, 5001, 5005};
Volume(60003) = {6003};
// 3-6 o'clock
Surface Loop(6004) = {5010, 5011, 4003, 3002, 5002, 5006};
Volume(60004) = {6004};
// 6-9 o'clock
Surface Loop(6005) = {5011, 5012, 4004, 3003, 5003, 5007};
Volume(60005) = {6005};
// 9-12 o'clock
Surface Loop(6006) = {5012, 5009, 4005, 3004, 5004, 5008};
Volume(60006) = {6006};

// Use hex meshing
Transfinite Line{:} = 20;
Transfinite Surface{:};
Recombine Surface{:};
Transfinite Volume{:};
