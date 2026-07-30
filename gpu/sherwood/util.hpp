// Common case information struct, to make device kernel arguments
// more readable.
typedef struct caseinfo {
    // Number of elements local to this MPI rank.
    dlong Nelements;
    // Offset between vector field components.
    dlong fieldOffset;
    // Total liquid volume.
    dfloat liquid_volume;
    // Total gas volume.
    dfloat gas_volume;
    // Global average density, for use in pressure gradient force term.
    dfloat rho_average;
    dfloat cumulative_c_sink;
    dfloat cumulative_dt;
} caseinfo_t;
