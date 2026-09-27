#!/usr/bin/env sh
set -eu

build_dir=${1:-build}
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
tester="$root/$build_dir/test/unit/tblite-tester"
tblite="$root/$build_dir/app/tblite"

if [ ! -x "$tester" ] || [ ! -x "$tblite" ]; then
  echo "Build artifacts not found in $build_dir" >&2
  exit 2
fi

for suite in \
  gfn1-hessian-validation \
  vibration \
  halogen \
  solvation-born \
  solvation-surface \
  solvation-cds
do
  "$tester" "$suite"
done

# Run the dedicated public C API smoke test.  The target is intentionally
# separate from the broad C API suite so this validation remains independent
# of optional ddX support.
capi=""
for candidate in \
  "$root/$build_dir/test/api/tblite-hessian-api-tester" \
  "$root/$build_dir/test/api/tblite-hessian-capi-tester"
do
  if [ -x "$candidate" ]; then
    capi=$candidate
    break
  fi
done
if [ -z "$capi" ]; then
  echo "GFN1 Hessian C API smoke tester not found in $build_dir" >&2
  exit 2
fi
"$capi"

tmp=${TMPDIR:-/tmp}/tblite-gfn1-hessian-$$
mkdir -p "$tmp"
trap 'rm -rf "$tmp"' EXIT HUP INT TERM

"$tblite" run --no-restart --method gfn1 --freq \
  --json "$tmp/gas.json" "$root/test/cli/10-water.xyz" >/dev/null
"$tblite" run --no-restart --method gfn1 --hess --alpb water \
  --json "$tmp/alpb.json" "$root/test/cli/10-water.xyz" >/dev/null

python3 - "$tmp/gas.json" "$tmp/alpb.json" <<'PY'
import json, math, sys
for fn in sys.argv[1:]:
    with open(fn, encoding="utf-8") as f:
        data = json.load(f)
    h = data.get("hessian")
    if not h or not all(math.isfinite(float(x)) for x in h):
        raise SystemExit(f"invalid Hessian in {fn}")
    if fn.endswith("gas.json"):
        freq = data.get("frequencies")
        if not freq or not all(math.isfinite(float(x)) for x in freq):
            raise SystemExit(f"invalid frequencies in {fn}")
PY

expect_hessian_failure() {
  name=$1
  shift
  if "$tblite" "$@" >/dev/null 2>"$tmp/$name.err"; then
    echo "$name Hessian unexpectedly succeeded" >&2
    exit 1
  fi
  if ! grep -qiE 'hessian|periodic|external|gfn1' "$tmp/$name.err"; then
    echo "$name failed without a Hessian-related diagnostic" >&2
    cat "$tmp/$name.err" >&2
    exit 1
  fi
}

# Unsupported combinations must fail explicitly rather than falling back to a
# numerical Hessian or silently dropping second-derivative terms.
expect_hessian_failure gfn2 run --no-restart --method gfn2 --hess \
  "$root/test/cli/10-water.xyz"
expect_hessian_failure external-field run --no-restart --method gfn1 --hess \
  --efield 0.1,0.0,0.0 "$root/test/cli/10-water.xyz"
expect_hessian_failure periodic run --no-restart --method gfn1 --hess \
  "$root/test/cli/08-gh131.gen"

echo "GFN1 analytic-Hessian release validation passed"
