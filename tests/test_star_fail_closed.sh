#!/usr/bin/env bash
# Exercise production STAR failure paths using tiny synthetic fixtures.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/scpolaseq-star-failure.XXXXXX")"
trap 'rm -rf "$TMP_ROOT"' EXIT

mkdir -p "$TMP_ROOT/bin" "$TMP_ROOT/star_index"
touch "$TMP_ROOT/star_index/SA" "$TMP_ROOT/star_index/Genome"
cat > "$TMP_ROOT/bin/STAR" <<'EOF'
#!/usr/bin/env bash
echo 'synthetic STAR failure' >&2
exit 42
EOF
chmod +x "$TMP_ROOT/bin/STAR"
export PATH="$TMP_ROOT/bin:$PATH"

run_expected_failure() {
    local name="$1" harness="$2" expected="$3"
    local work_dir="$TMP_ROOT/work_${name}"
    local log_file="$TMP_ROOT/${name}.log"
    local rc=0

    if nextflow run "$REPO_DIR/$harness" \
        -profile test \
        -work-dir "$work_dir" \
        --outdir "$TMP_ROOT/results_${name}" \
        --test_star_index "$TMP_ROOT/star_index" \
        -ansi-log false >"$log_file" 2>&1; then
        echo "FAIL: ${name} unexpectedly succeeded" >&2
        cat "$log_file" >&2
        return 1
    else
        rc=$?
    fi

    if [[ "$rc" -eq 0 ]] || ! find "$work_dir" -type f -name '.command.err' \
        -exec grep -F "$expected" {} + >/dev/null; then
        echo "FAIL: ${name} did not record expected failure '$expected' (exit ${rc})" >&2
        cat "$log_file" >&2
        find "$work_dir" -name '.command.err' -type f -print -exec cat {} \; >&2 || true
        return 1
    fi
    echo "PASS: ${name} failed closed as expected (exit ${rc})"
}

run_expected_failure alignment \
    tests/modules/alignment/test_starsolo_fail_closed.nf \
    'STAR alignment failed'
run_expected_failure reference \
    tests/modules/reference/test_prepare_reference_fail_closed.nf \
    'STAR index generation failed'
