# match-nek5000-202605

This case aims to match the Nek5000 run submitted in the ANS Winter Conference 2026 summary.
* Hard sink (c removed at interface and used to calculate interface transfer rate).
* Soft source (c added in liquid bulk to try to achieve quasi-steady-state).
* 32x64x32 mesh, order 5. Targets 50% Kolmogorov scale.
* dt=1e-4, CLSr=1e-3, TLSr=1e-2 intervals.
* Pressure residual tolerance and c SVV conditions are different from main case.

It has several changes from the original May 2026 Nek5000 conditions:
* The pressure preconditioner configuration does not work as-is in NekRS so that was commented out in the .par file. Not sure how to port that.
* The pressure `residualProj` was also commented out as it doesn't work as-is in NekRS.
* The soft source term psi*(1-c) was implemented as an explicit source term, because NekRS does not support implicit source terms. Nek5000 supported adding a psi explicit source and a -psi*c implicit source which is more stable.

The May 2026 Nek5000 conditions had some known issues:
* It assumed Fr as nu^2/(gL), not nu/sqrt(gL) but specified the Fr value (correctly) as nu/sqrt(gL). So it applies too much gravity. For parity this case maintains the same error.
* The 2x4x2 domain with 1 unit bubble diameter is too small and quickly leads to self-interaction and oscillation.
