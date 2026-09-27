#!/usr/bin/env python3
"""Build a pip-installable tblite wheel with Meson/meson-python.

The upstream Python binding is an independent Meson project that expects the
native tblite source as a Meson subproject.  This helper creates that layout in
a temporary directory without modifying the checkout and without requiring
symlink privileges on Windows.

The active Python interpreter is used as the wheel build interpreter.  Compiler
selection is controlled through the normal Meson environment variables, e.g.
FC=gfortran/CC=gcc or FC=ifx/CC=icx.
"""
from __future__ import annotations

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import platform
import configparser
import re


def _ignore_source(path: str, names: list[str]) -> set[str]:
    ignored: set[str] = set()
    fixed = {
        '.git', '.venv', 'venv', 'env', 'ENV', '__pycache__', '.pytest_cache',
        '.mypy_cache', '.ruff_cache', '.cache', 'dist', 'wheelhouse',
        '.venv-wheel', '.venv-install', '.wheel-build',
        'meson-build', 'builddir', '.fpm', 'fpm_build', 'htmlcov',
    }
    for name in names:
        if name in fixed:
            ignored.add(name)
            continue
        if name.startswith('.venv') or name.startswith('build-') or name.startswith('cmake-build-'):
            ignored.add(name)
            continue
        if name.endswith(('.pyc', '.pyo', '.whl')):
            ignored.add(name)
    return ignored


def _compiler_name(value: str | None) -> str:
    if not value:
        return ''
    return Path(value.strip('"\'')).name.lower()



def _select_lapack_backend(fc_name: str, requested: str | None) -> str:
    """Choose a Meson BLAS/LAPACK backend suitable for wheel builds.

    Windows ifx is reported by Meson as intel-llvm-cl.  dftd4 4.2.0 does
    not recognize that compiler id in its own automatic backend selection,
    so select MKL explicitly.  The wheel disables OpenMP, so the static MKL
    configuration uses mkl_sequential and avoids mkl_rt's runtime selection
    of the Intel threaded DLLs.  Windows GFortran builds are expected to use
    the MSYS2 OpenBLAS package.
    """
    if requested:
        return requested
    env_value = os.environ.get('TBLITE_LAPACK')
    if env_value:
        return env_value
    is_windows = os.name == 'nt' or platform.system().lower().startswith('win')
    if fc_name.startswith('ifx'):
        return 'mkl'
    if fc_name.startswith('gfortran') and is_windows:
        return 'openblas'
    return 'auto'


def _prepare_build_environment(fc_name: str) -> dict[str, str]:
    """Return an environment suitable for Meson native subprojects.

    Windows ifx does not reliably infer preprocessing from the .F90 suffix.
    Meson fallback dependencies such as mctc-lib and s-dftd3 contain .F90
    sources with C-preprocessor directives, so propagate /fpp through
    FFLAGS.  Meson recognizes FFLAGS as the Fortran compiler-flag environment
    variable and applies it to the full native build, including fallback
    subprojects.
    """
    env = os.environ.copy()
    is_windows = os.name == 'nt' or platform.system().lower().startswith('win')
    if is_windows and fc_name.startswith('ifx'):
        current = env.get('FFLAGS', '').strip()
        tokens = current.split() if current else []
        if not any(tok.lower() in {'/fpp', '-fpp'} for tok in tokens):
            env['FFLAGS'] = (current + ' /fpp').strip()
    return env



def _read_git_wrap(wrap_file: Path) -> tuple[str, str]:
    """Return (url, revision) from a Meson wrap-git file."""
    parser = configparser.ConfigParser()
    parser.read(wrap_file, encoding='utf-8')
    try:
        section = parser['wrap-git']
        return section['url'], section.get('revision', 'HEAD')
    except KeyError as exc:
        raise SystemExit(f'Invalid git wrap file: {wrap_file}') from exc


def _prepare_dftd4_static_for_windows_ifx(tblite_root: Path) -> None:
    """Prepare dftd4 for static linking with Windows ifx.

    dftd4 4.2.0 annotates many Fortran module procedures with Intel
    ``!DEC$ ATTRIBUTES DLLEXPORT`` directives.  When a client USEs those
    modules, Intel Fortran propagates the DLL property and emits references to
    ``__imp_*`` symbols.  That is correct for a DLL but incompatible with the
    static dftd4 library used by the wheel build.

    For the wheel only, fetch the exact revision from dftd4.wrap, remove the
    DLLEXPORT directives, and omit dftd4's standalone application/tests.  The
    library implementation itself is otherwise unchanged.
    """
    subprojects = tblite_root / 'subprojects'
    wrap_file = subprojects / 'dftd4.wrap'
    dftd4_dir = subprojects / 'dftd4'
    if dftd4_dir.exists():
        shutil.rmtree(dftd4_dir)

    url, revision = _read_git_wrap(wrap_file)
    print(f'[wheel] prefetching dftd4 {revision} for Windows ifx static build')
    subprocess.run(
        ['git', 'clone', '--depth', '1', '--branch', revision, url, str(dftd4_dir)],
        check=True,
    )

    removed = 0
    directive = re.compile(r'^\s*!DEC\$\s+ATTRIBUTES\s+DLLEXPORT\b.*$', re.IGNORECASE)
    for path in list(dftd4_dir.rglob('*.f90')) + list(dftd4_dir.rglob('*.F90')):
        text = path.read_text(encoding='utf-8')
        lines = text.splitlines()
        new_lines = [line for line in lines if not directive.match(line)]
        removed += len(lines) - len(new_lines)
        if len(new_lines) != len(lines):
            path.write_text('\n'.join(new_lines) + '\n', encoding='utf-8')

    meson_file = dftd4_dir / 'meson.build'
    meson_text = meson_file.read_text(encoding='utf-8')
    old_app = "# Add applications\nsubdir('app')"
    new_app = "# Add applications only for a standalone/shared build.\nif install\n  subdir('app')\nendif"
    old_test = "# add the testsuite\nsubdir('test')"
    new_test = "# Add the testsuite only for a standalone/shared build.\nif install\n  subdir('test')\nendif"
    if old_app not in meson_text or old_test not in meson_text:
        raise SystemExit('Unexpected dftd4 meson.build layout; static wheel patch was not applied')
    meson_text = meson_text.replace(old_app, new_app).replace(old_test, new_test)
    meson_file.write_text(meson_text, encoding='utf-8')

    config_file = dftd4_dir / 'config' / 'meson.build'
    config_text = config_file.read_text(encoding='utf-8')
    compiler_checks = "if fc_id == 'intel'"
    tbb_selection = "mkl_dep += cc.find_library('mkl_tbb_thread')"
    if config_text.count(compiler_checks) != 3 or tbb_selection not in config_text:
        raise SystemExit(
            'Unexpected dftd4 config/meson.build layout; static wheel patch was not applied'
        )
    config_text = config_text.replace(
        compiler_checks, "if fc_id == 'intel' or fc_id == 'intel-llvm-cl'"
    ).replace(tbb_selection, "mkl_dep += cc.find_library('mkl_sequential')")
    config_file.write_text(config_text, encoding='utf-8')
    shutil.rmtree(dftd4_dir / '.git', ignore_errors=True)
    print(
        f'[wheel] patched dftd4 static build: removed {removed} DLLEXPORT directives; '
        'app/tests disabled; enabled intel-llvm-cl; selected mkl_sequential'
    )


def _prepare_multicharge_static_for_windows_ifx(
    tblite_root: Path, project_root: Path
) -> None:
    """Prepare multicharge for the Windows ifx sequential-MKL wheel build.

    multicharge 0.5.0 only recognizes the classic ``intel`` Meson compiler
    identifier.  Current ifx is reported as ``intel-llvm-cl`` on Windows, so
    its MKL configuration otherwise rejects the explicit backend inherited
    from tblite.  The wheel has OpenMP disabled, therefore replace its stale
    TBB selection with mkl_sequential as well.
    """
    wrap_file = tblite_root / 'subprojects' / 'multicharge.wrap'
    multicharge_dir = project_root / 'subprojects' / 'multicharge'
    if multicharge_dir.exists():
        shutil.rmtree(multicharge_dir)

    url, revision = _read_git_wrap(wrap_file)
    print(f'[wheel] prefetching multicharge {revision} for Windows ifx static build')
    subprocess.run(
        ['git', 'clone', '--depth', '1', '--branch', revision, url, str(multicharge_dir)],
        check=True,
    )

    meson_file = multicharge_dir / 'config' / 'meson.build'
    meson_text = meson_file.read_text(encoding='utf-8')
    compiler_checks = "if fc_id == 'intel'"
    tbb_selection = "mkl_dep += cc.find_library('mkl_tbb_thread')"
    if meson_text.count(compiler_checks) != 3 or tbb_selection not in meson_text:
        raise SystemExit(
            'Unexpected multicharge config/meson.build layout; static wheel patch was not applied'
        )
    meson_text = meson_text.replace(
        compiler_checks, "if fc_id == 'intel' or fc_id == 'intel-llvm-cl'"
    ).replace(tbb_selection, "mkl_dep += cc.find_library('mkl_sequential')")
    meson_file.write_text(meson_text, encoding='utf-8')
    shutil.rmtree(multicharge_dir / '.git', ignore_errors=True)
    print('[wheel] patched multicharge static build: enabled intel-llvm-cl; selected mkl_sequential')

def _check_compilers() -> None:
    fc = os.environ.get('FC')
    cc = os.environ.get('CC')
    if not fc:
        raise SystemExit(
            'FC is not set. Set FC=gfortran or FC=ifx before building the wheel.'
        )
    fc_name = _compiler_name(fc)
    if not (fc_name.startswith('gfortran') or fc_name.startswith('ifx')):
        print(
            f'warning: FC={fc!r}; this repository is primarily tested with '
            'gfortran and Intel ifx.',
            file=sys.stderr,
        )
    if not cc:
        expected = 'gcc' if fc_name.startswith('gfortran') else 'icx'
        raise SystemExit(
            f'CC is not set. The Python CFFI extension also needs a C compiler; '
            f'use CC={expected} with FC={fc_name or fc}.'
        )
    for var in ('FC', 'CC'):
        exe = os.environ[var]
        if shutil.which(exe) is None and not Path(exe).exists():
            raise SystemExit(f'{var}={exe!r} was not found on PATH.')



def _existing_unique_dirs(paths: list[Path]) -> list[Path]:
    """Return existing directories without duplicates, preserving order."""
    seen: set[str] = set()
    result: list[Path] = []
    for path in paths:
        try:
            path = path.resolve()
        except OSError:
            continue
        if not path.is_dir():
            continue
        key = os.path.normcase(str(path))
        if key not in seen:
            seen.add(key)
            result.append(path)
    return result


def _windows_dll_search_dirs(fc_name: str) -> list[Path]:
    """Collect likely DLL directories for a Windows native wheel.

    delvewheel already searches PATH.  The extra paths here make locally built
    Intel oneAPI wheels robust even when the runtime/redist directories are not
    all present on PATH after the compiler process exits.
    """
    candidates: list[Path] = []

    # Existing PATH entries are useful for MSYS2/GFortran as well as oneAPI.
    for entry in os.environ.get('PATH', '').split(os.pathsep):
        if entry:
            candidates.append(Path(entry.strip('"')))

    mklroot = os.environ.get('MKLROOT')
    if mklroot:
        root = Path(mklroot)
        candidates.extend([
            root / 'bin',
            root / 'redist' / 'intel64',
            root / 'redist' / 'intel64_win',
        ])

    oneapi = os.environ.get('ONEAPI_ROOT')
    if oneapi:
        root = Path(oneapi)
        candidates.extend([
            root / 'compiler' / 'latest' / 'bin',
            root / 'compiler' / 'latest' / 'windows' / 'bin',
            root / 'compiler' / 'latest' / 'redist' / 'intel64_win' / 'compiler',
            root / 'compiler' / 'latest' / 'windows' / 'redist' / 'intel64_win' / 'compiler',
            root / 'mkl' / 'latest' / 'bin',
            root / 'mkl' / 'latest' / 'redist' / 'intel64',
        ])

    cmplr_root = os.environ.get('CMPLR_ROOT')
    if cmplr_root:
        root = Path(cmplr_root)
        candidates.extend([
            root / 'bin',
            root / 'windows' / 'bin',
            root / 'redist' / 'intel64_win' / 'compiler',
            root / 'windows' / 'redist' / 'intel64_win' / 'compiler',
        ])

    # Derive a compiler installation root directly from ifx/icx locations.
    if fc_name.startswith('ifx'):
        for exe_name in (os.environ.get('FC', 'ifx'), os.environ.get('CC', 'icx')):
            exe = shutil.which(exe_name)
            if not exe:
                continue
            exe_path = Path(exe).resolve()
            # Typical layout: <compiler>/<version>/bin/ifx.exe
            for parent in [exe_path.parent, *list(exe_path.parents)[:4]]:
                candidates.extend([
                    parent,
                    parent / 'redist' / 'intel64_win' / 'compiler',
                    parent / 'windows' / 'redist' / 'intel64_win' / 'compiler',
                ])

    return _existing_unique_dirs(candidates)


def _repair_windows_wheel(raw_wheel: Path, outdir: Path, fc_name: str) -> Path:
    """Vendor non-system DLL dependencies into a Windows wheel with delvewheel."""
    if shutil.which('delvewheel') is None:
        raise SystemExit(
            'delvewheel was not found. Reinstall the wheel build requirements with: '
            'python -m pip install -r tools/requirements-wheel-build.txt'
        )

    add_dirs = _windows_dll_search_dirs(fc_name)
    add_path = os.pathsep.join(str(path) for path in add_dirs)
    print(f'[wheel] repairing Windows DLL dependencies with delvewheel')
    if add_dirs:
        print('[wheel] delvewheel add-path directories:')
        for path in add_dirs:
            print('        ', path)

    show_cmd = ['delvewheel', 'show', str(raw_wheel)]
    repair_cmd = ['delvewheel', 'repair', str(raw_wheel), '-w', str(outdir)]
    if add_path:
        show_cmd += ['--add-path', add_path]
        repair_cmd += ['--add-path', add_path]

    # Diagnostic output is useful if a later import still fails.
    subprocess.run(show_cmd, check=True)
    existing = outdir / raw_wheel.name
    if existing.exists():
        existing.unlink()
    subprocess.run(repair_cmd, check=True)

    candidates = sorted(outdir.glob(raw_wheel.stem + '*.whl'), key=lambda p: p.stat().st_mtime)
    if not candidates:
        # delvewheel normally keeps the same wheel basename.
        candidates = sorted(outdir.glob('*.whl'), key=lambda p: p.stat().st_mtime)
    if not candidates:
        raise SystemExit('delvewheel completed but no repaired wheel was produced')
    repaired = candidates[-1]
    print('[wheel] repaired Windows wheel:', repaired)
    return repaired


def main() -> int:
    parser = argparse.ArgumentParser(
        description='Build a pip-installable wheel using Meson/meson-python.'
    )
    parser.add_argument(
        '--outdir', default='dist', help='wheel output directory (default: dist)'
    )
    parser.add_argument(
        '--keep-stage', action='store_true',
        help='keep the temporary staging directory for debugging',
    )
    parser.add_argument(
        '--lapack',
        choices=['auto', 'mkl', 'mkl-rt', 'openblas', 'netlib', 'custom'],
        default=None,
        help=(
            'BLAS/LAPACK backend forwarded to the native tblite Meson '
            'subproject. Default: mkl for ifx, openblas for Windows+gfortran, '
            'auto otherwise. '
            'TBLITE_LAPACK can also set the default.'
        ),
    )
    args = parser.parse_args()

    root = Path(__file__).resolve().parents[1]
    python_src = root / 'python'
    outdir = (root / args.outdir).resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    _check_compilers()
    try:
        import build  # noqa: F401
        import mesonpy  # noqa: F401
    except ImportError as exc:
        raise SystemExit(
            'Wheel build tools are missing. Install them first with: '
            'python -m pip install -r tools/requirements-wheel-build.txt'
        ) from exc
    for exe in ('meson', 'ninja'):
        if shutil.which(exe) is None:
            raise SystemExit(
                f'{exe} was not found on PATH. Install tools/requirements-wheel-build.txt first.'
            )

    # Use an ordinary temporary directory so a failed build can be kept for
    # inspection.  On success it is removed unless --keep-stage was used.
    stage = Path(tempfile.mkdtemp(prefix='tblite-wheel-'))

    project = stage / 'python'
    print(f'[wheel] staging project in {stage}')
    shutil.copytree(
        python_src,
        project,
        ignore=shutil.ignore_patterns(
            'build', 'dist', '*.egg-info', '__pycache__', '.pytest_cache',
            'subprojects',
        ),
    )
    subproject = project / 'subprojects' / 'tblite'
    subproject.parent.mkdir(parents=True, exist_ok=True)
    shutil.copytree(root, subproject, ignore=_ignore_source)

    fc_name = _compiler_name(os.environ.get('FC'))
    lapack = _select_lapack_backend(fc_name, args.lapack)
    is_windows = os.name == 'nt' or platform.system().lower().startswith('win')
    if is_windows and fc_name.startswith('ifx'):
        if shutil.which('git') is None:
            raise SystemExit(
                'git was not found on PATH. Windows ifx wheel builds prefetch '
                'the pinned dftd4 source so it can be prepared for static linking.'
            )
        _prepare_dftd4_static_for_windows_ifx(subproject)
        _prepare_multicharge_static_for_windows_ifx(subproject, project)
    raw_outdir = stage / 'raw-wheel' if is_windows else outdir
    raw_outdir.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable, '-m', 'build', '--wheel', '--no-isolation',
        '-Csetup-args=-Dtblite:lapack=' + lapack,
        '--outdir', str(raw_outdir), str(project),
    ]
    build_env = _prepare_build_environment(fc_name)
    print('[wheel] FC =', os.environ.get('FC'))
    print('[wheel] CC =', os.environ.get('CC'))
    print('[wheel] FFLAGS =', build_env.get('FFLAGS', ''))
    print('[wheel] LAPACK backend =', lapack)
    print('[wheel] running:', ' '.join(cmd))
    try:
        subprocess.run(cmd, check=True, env=build_env)
    except subprocess.CalledProcessError:
        # Keep the complete Meson project and logs on failure.
        print('[wheel] build failed; staging directory kept at:', stage, file=sys.stderr)
        logs = list(project.glob('.mesonpy-*/meson-logs/meson-log.txt'))
        if logs:
            print('[wheel] Meson log:', logs[0], file=sys.stderr)
        if os.name == 'nt' and fc_name.startswith('ifx'):
            print(
                '[wheel] Windows ifx requires an initialized oneAPI environment '
                'and Fortran preprocessing. This helper sets FFLAGS=/fpp '
                'automatically; also verify MKLROOT after running setvars.bat.',
                file=sys.stderr,
            )
        if os.name == 'nt' and fc_name.startswith('gfortran'):
            print(
                '[wheel] Windows gfortran defaults to OpenBLAS. Install the '
                'matching MSYS2 OpenBLAS package or override --lapack.',
                file=sys.stderr,
            )
        raise

    raw_wheels = sorted(raw_outdir.glob('*.whl'), key=lambda p: p.stat().st_mtime)
    if not raw_wheels:
        raise SystemExit('build completed but no wheel was produced')
    raw_wheel = raw_wheels[-1]

    if is_windows:
        final_wheel = _repair_windows_wheel(raw_wheel, outdir, fc_name)
    else:
        final_wheel = raw_wheel

    print('[wheel] created:', final_wheel)
    if args.keep_stage:
        print('[wheel] staging directory kept at:', stage)
    else:
        shutil.rmtree(stage, ignore_errors=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
