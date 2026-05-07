// Common case information struct, to make device kernel arguments
// more readable.
typedef struct caseinfo {
    // Number of elements local to this CPU node.
    dlong Nelements;
    // Offset between vector field components.
    dlong fieldOffset;
    // Global average density, for use in pressure gradient force term.
    dfloat rho_average;
} caseinfo_t;
