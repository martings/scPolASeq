#!/usr/bin/env bash
# Only Nextflow workflow preflight is exercised; no downstream task may start.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/scpolaseq-downstream-preflight.XXXXXX")"
trap 'rm -rf -- "$TEST_DIR"' EXIT

run_rejected() {
    local name="$1" python_exec="$2" output="$3" expected="$4"
    if nextflow run "$REPO_DIR/workflows/pbmc1k_downstream.nf" \
        --upstream_results "$TEST_DIR" \
        --downstream_package_src "$TEST_DIR" \
        --downstream_outdir "$output" \
        --downstream_python "$python_exec" \
        -work-dir "$TEST_DIR/work-$name" -ansi-log false \
        > "$TEST_DIR/$name.log" 2>&1; then
        echo "FAIL: $name was accepted" >&2
        return 1
    fi
    grep -F "$expected" "$TEST_DIR/$name.log" >/dev/null || {
        cat "$TEST_DIR/$name.log" >&2
        return 1
    }
    if grep -F 'Submitted process' "$TEST_DIR/$name.log" >/dev/null; then
        echo "FAIL: $name submitted a process" >&2
        return 1
    fi
}

unsafe_target="$TEST_DIR/unsafe-target"
run_rejected unsafe_python 'python3;id' "$unsafe_target" \
    'Invalid --downstream_python'
test ! -e "$unsafe_target"

existing_target="$TEST_DIR/existing-target"
mkdir "$existing_target"
run_rejected existing_output python3 "$existing_target" \
    'Downstream output already exists'
test -d "$existing_target"
echo 'PASS: downstream preflight rejects shell-shaped Python and existing output'
