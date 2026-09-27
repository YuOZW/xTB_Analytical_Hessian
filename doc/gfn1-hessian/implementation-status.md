# GFN1-xTB analytic Hessian implementation status

The analytic Cartesian Hessian implementation is complete for the intended
molecular GFN1-xTB scope and has reached the upstream-integration stage.

## Supported

- non-periodic GFN1-xTB
- closed-shell and collinear open-shell calculations
- zero and finite electronic temperature
- SCC2/SCC3 response
- D3(BJ), coordination-dependent H0, repulsion, and GFN1 halogen correction
- ALPB/GBSA solvation and supported CDS contributions
- Fortran, C, Python, QCSchema, and CLI access
- harmonic frequencies, projected normal modes, and imaginary-mode counting

## Intentionally unsupported

- GFN2-xTB analytic Hessians
- periodic Hessians / phonons
- external electric fields
- ddX solvation Hessians

Unsupported combinations fail explicitly.

## Validation

`tools/validate-gfn1-hessian.sh <build-dir>` is the release validation entry
point.  It checks broad finite-difference agreement against the analytic
gradient, Hessian symmetry, translational invariance, rotational covariance,
halogen and solvation terms, vibrational analysis, public C API access,
gas-phase and ALPB CLI paths, and explicit rejection of GFN2, periodic and
external-field Hessians.

The broad validation covers neutral, charged, open-shell, finite-temperature,
degenerate-frontier, diatomic, main-group, and chemically diverse systems.
Exact occupied/virtual frontier degeneracy at zero electronic temperature is
rejected when the state-specific response is not differentiable; finite
electronic temperature can be used for an appropriate state-averaged response.

## Integration state

Phase 19 rebuilt the deliverable from the pristine upstream source snapshot,
removed generated artifacts, registered a standalone C API Hessian smoke test
in CMake and Meson, and generated a clean upstream patch.  The implementation
is therefore ready for commit splitting / PR preparation rather than further
Hessian feature development.
