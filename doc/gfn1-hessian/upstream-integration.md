# GFN1-xTB analytic Hessian: upstream integration notes

This tree is a clean integration snapshot based on tblite main commit
`c701d0201b749975bd70ff56225abb074a4f6372` (source artifact dated
2026-09-18).  It contains only source, tests, documentation and validation
assets required for the GFN1-xTB analytic-Hessian implementation; build and
runtime-generated files have been removed.

## Supported scope

- GFN1-xTB, non-periodic molecular calculations
- closed- and collinear open-shell calculations
- zero and finite electronic temperature
- D3(BJ), SCC2/SCC3, coordination-dependent H0, repulsion and halogen correction
- ALPB/GBSA and CDS solvation
- Cartesian Hessian, harmonic frequencies and normal modes
- Fortran, C and Python APIs, CLI and QCSchema Hessian driver

Explicitly unsupported for analytic Hessians:

- GFN2-xTB and other Hamiltonians
- periodic systems
- external electric fields
- ddX solvation

Unsupported combinations fail explicitly; they never silently fall back to a
numerical Hessian.

## Build-system integration

The new Fortran sources are registered in both CMake and Meson manifests:

- `src/tblite/ncoord.f90`
- `src/tblite/solvation/hyperdual.f90`
- `src/tblite/xtb/response.f90`
- `src/tblite/xtb/hessian.f90`
- `src/tblite/xtb/vibration.f90`

A standalone public-C-API smoke target is also registered in both CMake and
Meson.  `fpm.toml` is unchanged because the existing fpm source discovery
covers the new files without an explicit source list.

## Clean validation

Configure a fresh CMake build with tests and the C API enabled, then run:

```sh
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DTBLITE_WITH_TESTS=ON \
  -DTBLITE_WITH_API=ON \
  -DTBLITE_WITH_DDX=OFF
cmake --build build
./tools/validate-gfn1-hessian.sh build
```

The release validation checks:

- analytic-vs-numerical-gradient Hessian validation,
- vibrational analysis,
- halogen and ALPB/GBSA/CDS second derivatives,
- the standalone public C Hessian API,
- gas-phase and ALPB CLI Hessian paths,
- explicit rejection of GFN2, periodic and external-field Hessians.

For Phase 19 the fresh Release build completed successfully.  The dedicated C
API test, vibration, halogen, Born/ALPB, surface and CDS suites, positive CLI
paths, and all three unsupported-mode rejection checks passed.  A rerun of the
broad Phase-13 validation passed its diatomic, main-group, charged, open-shell
and frontier-degeneracy groups; the final MB16-43 numerical-reference group is
long-running and exceeded this execution environment's per-command time cap.
That test and the underlying Hessian code are unchanged from Phase 18, where
the complete broad validation was already run to completion.

## Patch

`GFN1_HESSIAN_UPSTREAM.patch` is generated against the pristine upstream
snapshot and can be checked/applied with:

```sh
git apply --check GFN1_HESSIAN_UPSTREAM.patch
git apply GFN1_HESSIAN_UPSTREAM.patch
```

The patch intentionally does not contain build products or runtime-generated
files.
