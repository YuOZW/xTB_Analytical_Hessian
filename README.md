# tblite 0.7.0 GFN1-xTB Analytic Hessian

This repository is a modified fork of upstream [`tblite` 0.7.0](https://github.com/tblite/tblite) that adds a fully analytic Cartesian Hessian and harmonic vibrational analysis for molecular GFN1-xTB calculations.

The implementation differentiates the analytic GFN1-xTB gradient and solves the first-order self-consistent electronic response equations. Production Hessians therefore do **not** use finite-difference derivatives.

This repository is **not an official upstream tblite release** and is not endorsed by the tblite or xTB developers.

## Features

- Analytic Cartesian Hessian for non-periodic GFN1-xTB
- Closed-shell and collinear open-shell calculations
- Zero and finite electronic temperature
- Analytic response of SCC2 and SCC3 terms
- Coordination-dependent H0, repulsion, and D3(BJ) contributions
- GFN1 halogen-correction Hessian
- ALPB and GBSA solvation, including supported CDS contributions
- Harmonic vibrational frequencies and Cartesian normal modes
- Translational and rotational projection for vibrational analysis
- Imaginary-frequency counting for minimum / transition-state characterization
- Fortran, C, Python, CLI, and QCSchema interfaces
- Validation against central differences of the analytic gradient

See [`doc/properties/hessian.rst`](doc/properties/hessian.rst) and [`doc/gfn1-hessian/implementation-status.md`](doc/gfn1-hessian/implementation-status.md) for implementation details and the current validation scope.

## Important Notes

This is research software under active development. **Numerical accuracy, completeness, reproducibility, and suitability for scientific or production use are not guaranteed.** Results should be independently validated before they are used for scientific conclusions.

The analytic Hessian currently supports **GFN1-xTB only**. The following are intentionally unsupported:

- GFN2-xTB analytic Hessians
- Periodic systems and phonons
- External electric fields
- ddX solvation Hessians

Unsupported combinations fail explicitly rather than silently omitting a second-derivative contribution.

At zero electronic temperature, an exact occupied/virtual frontier degeneracy can make the state-specific electronic response non-unique. In such cases the Hessian routine reports an error. A small finite electronic temperature can be used when a state-averaged response is scientifically appropriate.

### Development disclosure

**OpenAI ChatGPT was used during the implementation of the GFN1-xTB analytic Hessian extension**, including implementation planning, code generation, refactoring, debugging, test generation, validation design, and documentation. The resulting code remains subject to human review and independent scientific validation. Use of ChatGPT does **not** imply that OpenAI has reviewed, validated, endorsed, or supports this software or its scientific results.

### Upstream relationship

This work is derived from [`tblite`](https://github.com/tblite/tblite). The analytic-Hessian implementation was developed against upstream commit:

```text
c701d0201b749975bd70ff56225abb074a4f6372
```

The original tblite project and its developers are not responsible for modifications introduced by this repository.

## Repository Layout

The upstream tblite directory layout is retained so Meson, CMake, Fortran
Package Manager, and the Python packaging configuration continue to work
without path-specific patches.

| Path | Purpose |
| --- | --- |
| `src/` | Fortran implementation, including the GFN1-xTB analytic-Hessian code |
| `include/` | Public C API headers |
| `python/` | Python/CFFI interface and wheel packaging project |
| `app/` | Native command-line application |
| `test/` | Native unit, API, CLI, and validation tests |
| `doc/` | Sphinx documentation and analytic-Hessian implementation notes |
| `tools/` | Build, validation, and wheel smoke-test helpers |
| `subprojects/` | Pinned Meson fallback dependency definitions |
| `.github/` | Continuous-integration and regression workflows |

Generated build directories, Python caches, and wheel artifacts under `dist/`
are deliberately ignored by Git. Build a release wheel locally or attach it to
a GitHub Release; do not commit it to the source repository.

## Requirements

For a native build:

- A Fortran 2008 compiler (GFortran or Intel oneAPI `ifx`)
- CMake >= 3.14 or Meson >= 0.57.2
- Ninja or another supported build backend
- BLAS / LAPACK

Optional components:

- A C compiler for the C API tests
- Python and CFFI for the Python interface
- A C++ compiler for upstream ddX functionality; ddX Hessians are not supported by this extension

## Installation

### Build a pip-installable wheel with Meson

The recommended Python installation flow is:

1. prepare a **build environment** with a Fortran/C compiler and Meson tools;
2. build a platform-specific wheel with Meson through `meson-python`;
3. create a separate **runtime Python environment**;
4. install the generated `.whl` with `pip`.

Meson performs the native Fortran/C build, while `meson-python` packages the Meson install result as a standard Python wheel. The wheel contains the Python `tblite` module and the statically linked native tblite library used by the CFFI extension. The supported GFN1-xTB Hessian wheel intentionally disables ddX, MPI, ScaLAPACK, HDF5, TREXIO, and OpenMP to reduce external runtime dependencies. ALPB/GBSA/CDS solvation used by the analytic Hessian remains available. On Windows, `tools/build-wheel.py` additionally runs `delvewheel` after the Meson build so non-system DLL dependencies of `_libtblite.pyd` (for example Intel Fortran/oneMKL runtime DLLs or GNU runtime DLLs) are vendored into the final wheel.

The wheel installs the **Python interface**. It does not install the native `tblite` command-line executable on `PATH`; build/install the native Meson target separately if the CLI is required.

The wheel is **platform- and Python-version-specific**. Build the wheel on the same OS/architecture family on which it will be installed. A wheel built on Linux cannot be installed on Windows, and vice versa.

The first build may download Meson fallback dependencies referenced by the files under `subprojects/`. An internet connection is therefore normally required unless those dependencies are already available in the Meson package cache.

#### 1. Create the wheel-build Python environment

Linux:

```bash
python3 -m venv .venv-wheel
source .venv-wheel/bin/activate
python -m pip install --upgrade pip
python -m pip install -r tools/requirements-wheel-build.txt
```

Windows PowerShell:

```powershell
py -m venv .venv-wheel
.\.venv-wheel\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r tools\requirements-wheel-build.txt
```

#### 2A. Linux with GFortran

A C compiler and BLAS/LAPACK development libraries are also required. For Debian/Ubuntu, for example:

```bash
sudo apt-get install gfortran gcc libblas-dev liblapack-dev
export FC=gfortran
export CC=gcc
python tools/build-wheel.py

# Optional explicit backend selection:
# python tools/build-wheel.py --lapack auto
```

The wheel is written to `dist/`.

#### 2B. Linux with Intel oneAPI `ifx`

Use the LLVM-based Intel compiler **`ifx`**, not `ifort`. Initialize oneAPI first:

```bash
source /opt/intel/oneapi/setvars.sh
export FC=ifx
export CC=icx
python tools/build-wheel.py

# On Linux+ifx the helper selects MKL automatically.
# To override it explicitly:
# python tools/build-wheel.py --lapack mkl
```

With `ifx`, Meson automatically selects Intel MKL for BLAS/LAPACK when the oneAPI environment is initialized.

#### 2C. Windows with GFortran (MSYS2 / MinGW-w64)

The recommended GFortran setup on Windows is **MSYS2 UCRT64** (or MinGW64) with a 64-bit Python interpreter. Install MSYS2 first, then open the **MSYS2 UCRT64** terminal and install the compiler toolchain:

```bash
pacman -Syu
pacman -S --needed mingw-w64-ucrt-x86_64-gcc \
  mingw-w64-ucrt-x86_64-gcc-fortran \
  mingw-w64-ucrt-x86_64-openblas
```

After installation, confirm that the compilers are available:

```bash
gfortran --version
gcc --version
```

For wheel builds, it is convenient to use the normal Windows Python installation from PowerShell and prepend the MSYS2 compiler directory to `PATH`. For the default UCRT64 installation:

```powershell
$env:Path = "C:\msys64\ucrt64\bin;" + $env:Path
$env:FC = "gfortran"
$env:CC = "gcc"

python --version
gfortran --version
gcc --version
python tools\build-wheel.py
```

For Windows+GFortran, `build-wheel.py` selects `openblas` automatically. You can also request it explicitly with:

```powershell
python tools\build-wheel.py --lapack openblas
```

If MinGW64 is used instead of UCRT64, use `C:\msys64\mingw64\bin` and the corresponding `mingw-w64-x86_64-*` packages. The Python interpreter, compiler toolchain, and target architecture must match; normally all components should be **64-bit**.

If Meson cannot locate BLAS/LAPACK automatically, install the matching MSYS2 BLAS/LAPACK packages or provide an appropriate Meson native file / linker search path. Do not mix 32-bit and 64-bit libraries.

#### 2D. Windows with Intel oneAPI `ifx` (recommended)

For Intel builds, use the LLVM-based **Intel Fortran Compiler `ifx`** rather than the classic `ifort`. Install:

1. **Intel oneAPI Base Toolkit** (provides `icx` and MKL),
2. **Intel oneAPI HPC Toolkit** (provides `ifx`), and
3. **Microsoft Visual Studio 2022 Build Tools** or Visual Studio with the **Desktop development with C++** workload.

The simplest option is to launch an **Intel oneAPI command prompt for Intel 64** from the Start menu. This initializes Visual Studio, `ifx`, `icx`, and MKL automatically. Confirm the environment before building:

```bat
ifx --version
icx --version
where ifx
where icx
```

Then create/activate the wheel-build Python environment and run:

```bat
set FC=ifx
set CC=icx
python tools\build-wheel.py

rem Equivalent explicit backend selection:
rem python tools\build-wheel.py --lapack mkl
```

When using a normal `cmd.exe`, initialize oneAPI first (the exact path can vary with the installed oneAPI version):

```bat
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
set FC=ifx
set CC=icx
python tools\build-wheel.py
```

In PowerShell, run the oneAPI environment through `cmd.exe` or start PowerShell from an initialized oneAPI command prompt. After initialization, set:

```powershell
$env:FC = "ifx"
$env:CC = "icx"
ifx --version
icx --version
$env:MKLROOT
python tools\build-wheel.py

# Equivalent explicit backend selection:
# python tools\build-wheel.py --lapack mkl
```

`CXX` is not required for the supported Hessian wheel because ddX is disabled in this wheel configuration. With oneAPI initialized, Windows+`ifx` wheel builds explicitly select static MKL (`mkl`). Because the wheel disables OpenMP, this links the `mkl_sequential` layer rather than `mkl_rt`'s dynamically selected Intel threaded layer. This avoids both the dftd4 4.2.0 Meson issue where Windows `ifx` is reported as `intel-llvm-cl` and the runtime failure caused by a missing `mkl_intel_thread.2.dll`. The wheel helper still runs `delvewheel` to bundle any remaining non-system compiler-runtime DLLs.

Windows `ifx` also requires explicit Fortran preprocessing for `.F90` sources used by tblite and several Meson fallback dependencies (`mctc-lib`, `s-dftd3`, etc.). `tools/build-wheel.py` now adds `/fpp` to `FFLAGS` automatically for Windows+`ifx`. `FFLAGS` is the Meson-recognized Fortran compiler flag variable and is inherited by the fallback subprojects (`mctc-lib`, `s-dftd3`, etc.). You should see a line similar to:

```text
[wheel] FFLAGS = /fpp
```

Do not remove this flag: without it, directives such as `#include`, `#ifdef`, and `#define` are parsed as invalid Fortran and symbols such as `WITH_JSON` or `TBLITE_HAS_DDX` remain undefined.

If a Windows `ifx` build reports `C shared or static library 'blas' not found`, update to the current `tools/build-wheel.py`, initialize oneAPI with `setvars.bat`, verify that `MKLROOT` is defined, and retry:

```powershell
$env:FC = "ifx"
$env:CC = "icx"
$env:MKLROOT
python tools\build-wheel.py --lapack mkl
```

The helper now also keeps its temporary staging directory after a failed build and prints the path to `meson-log.txt`, which makes dependency-detection failures easier to inspect.

If the build instead fails with messages such as `Bad # preprocessor line`, `WITH_JSON has no type`, `DFTD3_HAS_MPI has already been assigned a data type`, or undefined `TBLITE_HAS_*` constants, the Fortran preprocessor was not enabled. With the current helper this is automatic. Confirm that the build header contains:

```text
[wheel] FC = ifx
[wheel] CC = icx
[wheel] FFLAGS = /fpp
[wheel] LAPACK backend = mkl
```

If the build reaches the link stage and fails while creating `dftd4.exe` with errors such as `LNK2019: unresolved external symbol __imp_DFTD4_*`, the problem is different from MKL or preprocessing. dftd4 4.2.0 contains Intel `!DEC$ ATTRIBUTES DLLEXPORT` directives. When Windows `ifx` consumes those module files, it treats the exported procedures as DLL imports in downstream code, which is incompatible with the static dftd4 library used inside the wheel. The current `tools/build-wheel.py` handles this automatically for **Windows + ifx** by:

1. prefetching the exact dftd4 revision pinned by `subprojects/dftd4.wrap`,
2. removing `DLLEXPORT` directives only in the temporary wheel staging copy,
3. disabling the standalone dftd4 CLI and test targets for that static subproject, and
4. linking the patched static dftd4 library into tblite.

The upstream dftd4 source in your repository is not modified. At the start of the build you should see messages similar to:

```text
[wheel] prefetching dftd4 v4.2.0 for Windows ifx static build
[wheel] patched dftd4 static build: removed <N> DLLEXPORT directives; app/tests disabled
```

`git` must therefore be available on `PATH` for Windows+ifx wheel builds. If you previously saw `__imp_DFTD4_*` linker errors, use the current `tools/build-wheel.py` and rebuild from a fresh invocation; the helper always creates a new staging directory.

If the final Python-extension link fails with `LNK1104: cannot open file 'python3.lib'` while the same command already contains a version-specific library such as `python313.lib`, this comes from CFFI enabling CPython's Limited API on Windows. Conda environments commonly provide `python313.lib` (or the corresponding version-specific import library) but not `python3.lib`. This project builds normal CPython-version-specific wheels such as `cp313-win_amd64`, not `abi3` wheels, so the CFFI compilation explicitly defines `_CFFI_NO_LIMITED_API`. The resulting extension links against the version-specific Python import library only. No separate `python3.lib` installation or manual copy is required.

For a **native** Windows Meson build performed without `tools/build-wheel.py`, set the flag yourself before `meson setup` so it also reaches fallback subprojects:

```bat
set FC=ifx
set CC=icx
set FFLAGS=/fpp
meson setup build-ifx --buildtype=release -Dlapack=mkl
meson compile -C build-ifx
```

PowerShell equivalent:

```powershell
$env:FC = "ifx"
$env:CC = "icx"
$env:FFLAGS = "/fpp"
meson setup build-ifx --buildtype=release -Dlapack=mkl
meson compile -C build-ifx
```



##### Windows import error: `_libtblite` DLL load failed

If `pip install` succeeds but

```python
from tblite.interface import Calculator
```

fails with

```text
ImportError: DLL load failed while importing _libtblite
```

then `_libtblite.pyd` was found but one of its native DLL dependencies was not. This is especially common for an unrepaired Windows `ifx` wheel because the extension can depend on Intel Fortran runtime libraries and the oneMKL runtime.

The current `tools/build-wheel.py` handles this automatically. Make sure the wheel-build environment includes `delvewheel`:

```bat
python -m pip install -r tools\requirements-wheel-build.txt
```

Then rebuild from an initialized oneAPI environment:

```bat
call "C:\Program Files (x86)\Intel\oneAPI\setvars.bat"
set FC=ifx
set CC=icx
python tools\build-wheel.py
```

Near the end of the build you should see messages similar to:

```text
[wheel] repairing Windows DLL dependencies with delvewheel
[wheel] repaired Windows wheel: ...\dist\tblite-...-win_amd64.whl
```

Install that repaired wheel into the runtime environment:

```bat
python -m pip uninstall -y tblite
python -m pip install --force-reinstall dist\tblite-0.7.0-*.whl
python -c "from tblite.interface import Calculator; print('tblite import OK')"
```

To inspect an already-built Windows wheel manually, use:

```bat
delvewheel show dist\tblite-0.7.0-*.whl
```

`delvewheel show` lists external DLL dependencies, while `delvewheel repair` copies those DLLs into the wheel and patches the wheel's DLL-loading path. If you are troubleshooting a custom oneAPI installation, run `setvars.bat` first so the Intel compiler and MKL runtime directories are visible on `PATH`.


#### 3. Install the wheel in a separate Python environment

After the wheel has been generated, the build environment is no longer required. Create a clean runtime environment and install the wheel with `pip`.

Linux:

```bash
deactivate  # if the build venv is active
python3 -m venv .venv-install
source .venv-install/bin/activate
python -m pip install --upgrade pip
python -m pip install dist/tblite-0.7.0-*.whl
python tools/smoke-test-wheel.py
```

Windows PowerShell:

```powershell
deactivate  # if the build venv is active
py -m venv .venv-install
.\.venv-install\Scripts\Activate.ps1
python -m pip install --upgrade pip
# Replace <wheel-file> with the wheel created in dist\.
python -m pip install .\dist\<wheel-file>.whl
python tools\smoke-test-wheel.py
```

The runtime environment does not need Meson, Ninja, or a Fortran compiler. On Windows, `tools/build-wheel.py` automatically repairs the wheel with `delvewheel`, which copies the detected non-system DLL dependencies into the wheel and patches the package so they are discoverable at import time. Therefore the **final wheel written to `dist/`** should be installed, not the raw intermediate wheel from the temporary staging directory. On Linux, compiler/runtime portability is still the packager's responsibility; inspect and repair/bundle shared-library dependencies as appropriate (for example with `auditwheel` for manylinux distribution).

> **Note:** `python tools/build-wheel.py` is a convenience wrapper. It stages the existing `python/` Meson project together with this repository as a local Meson subproject and then runs `python -m build --wheel --no-isolation`. No symlink privileges are required, which keeps the same workflow usable on Windows and Linux. The helper chooses `mkl` for ifx, `openblas` for Windows+GFortran, and `auto` otherwise. Override this with `--lapack <backend>` or the `TBLITE_LAPACK` environment variable.

### Native build with Meson

For a native library/CLI build without creating a Python wheel:

```bash
meson setup build --buildtype=release
meson compile -C build
meson test -C build --print-errorlogs
```

For Intel oneAPI:

```bash
source /opt/intel/oneapi/setvars.sh
export FC=ifx
export CC=icx
meson setup build-ifx --buildtype=release
meson compile -C build-ifx
```

### Native build with CMake

```bash
cmake -S . -B build -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DTBLITE_WITH_TESTS=ON \
  -DTBLITE_WITH_API=ON
cmake --build build
```

Run the GFN1-xTB Hessian validation suite with:

```bash
./tools/validate-gfn1-hessian.sh build
```

## Quick Start

### Native command line

After a native Meson/CMake installation, evaluate the analytic Cartesian Hessian:

```bash
tblite run --method gfn1 --hess water.xyz
```

Evaluate the Hessian and harmonic vibrational frequencies:

```bash
tblite run --method gfn1 --freq water.xyz
```

Run a solvated calculation with ALPB:

```bash
tblite run --method gfn1 --freq --alpb water water.xyz
```

The Cartesian Hessian is reported in Hartree/Bohr². Harmonic frequencies are reported in cm⁻¹, with imaginary modes represented by negative frequencies.

### Python

```python
import numpy as np
from tblite.interface import Calculator

numbers = np.array([8, 1, 1])
positions = np.array([
    [ 0.0000000000, 0.0000000000, 0.0],
    [ 1.4323367300, 1.1071526600, 0.0],
    [-1.4323367300, 1.1071526600, 0.0],
])

calc = Calculator("GFN1-xTB", numbers, positions)
calc.set("temperature", 0.0)

res = calc.hessian()

hessian = res.get("hessian")
frequencies = res.get("frequencies")
modes = res.get("normal-modes")
n_imaginary = res.get("n-imaginary")
```

The Hessian has shape `(3N, 3N)`. The vibrational analysis automatically projects five rigid-body modes for linear molecules and six for nonlinear molecules.

## Validation

The release validation suite checks:

- analytic Hessian vs. central differences of the analytic gradient
- Hessian symmetry
- translational invariance
- rotational covariance
- coupled-response residuals
- closed-shell and open-shell systems
- finite electronic temperature
- charged and degenerate-frontier systems
- halogen correction
- ALPB / GBSA / CDS solvation
- harmonic vibrational analysis
- public C API access
- CLI Hessian and frequency paths
- explicit rejection of unsupported Hessian modes

The dedicated entry point is:

```bash
./tools/validate-gfn1-hessian.sh build
```

See [`doc/gfn1-hessian/upstream-integration.md`](doc/gfn1-hessian/upstream-integration.md) for integration and validation details.

## References

This extension builds on the GFN1-xTB method and the upstream tblite implementation. Relevant references include:

1. S. Grimme, C. Bannwarth, and P. Shushkov, **A Robust and Accurate Tight-Binding Quantum Chemical Method for Structures, Vibrational Frequencies, and Noncovalent Interactions of Large Molecular Systems Parametrized for All spd-Block Elements (Z = 1–86)**, *J. Chem. Theory Comput.* **13** (2017) 1989–2009. https://doi.org/10.1021/acs.jctc.7b00118

2. C. Bannwarth, E. Caldeweyher, S. Ehlert, A. Hansen, P. Pracht, J. Seibert, S. Spicher, and S. Grimme, **Extended tight-binding quantum chemistry methods**, *WIREs Comput. Mol. Sci.* **11** (2021) e1493. https://doi.org/10.1002/wcms.1493

If this implementation is used in scientific work, please cite the original GFN1-xTB methodological literature and the upstream tblite project as appropriate. A dedicated publication for this analytic-Hessian implementation can be added here if one becomes available.

## License

The **original GFN1-xTB analytic-Hessian contributions introduced by this repository are made available under the MIT License**. See [`LICENSE`](LICENSE).

This repository also contains source code derived from upstream `tblite`. Upstream tblite code remains subject to its original **GNU Lesser General Public License v3.0 or later (LGPL-3.0-or-later)** terms. See [`COPYING`](COPYING) and [`COPYING.LESSER`](COPYING.LESSER).

Because this repository is a derivative work containing LGPL-licensed upstream code, redistribution of the combined source must continue to satisfy the applicable LGPL requirements. The MIT license does not replace or remove the upstream tblite license obligations.
