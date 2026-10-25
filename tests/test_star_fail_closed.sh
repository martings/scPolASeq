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
case "${STAR_MODE:-fail}" in
    fail)
        echo 'synthetic STAR failure' >&2
        exit 42
        ;;
    matrix_missing)
        printf 'synthetic BAM\n' > fail_closed.Aligned.sortedByCoord.out.bam
        ;;
    bai_missing)
        printf 'synthetic BAM\n' > fail_closed.Aligned.sortedByCoord.out.bam
        mkdir -p fail_closed.Solo.out/Gene/filtered
        printf 'synthetic matrix\n' > fail_closed.Solo.out/Gene/filtered/matrix.mtx
        ;;
    success) ;;
    *) exit 64 ;;
esac
exit 0
EOF
cat > "$TMP_ROOT/bin/samtools" <<'EOF'
#!/usr/bin/env bash
# Simulate samtools success without creating its requested BAI output.
exit 0
EOF
chmod +x "$TMP_ROOT/bin/STAR"
chmod +x "$TMP_ROOT/bin/samtools"
export PATH="$TMP_ROOT/bin:$PATH"

run_expected_failure() {
    local name="$1" mode="$2" harness="$3" expected="$4"
    local work_dir="$TMP_ROOT/work_${name}"
    local log_file="$TMP_ROOT/${name}.log"
    local rc=0

    if STAR_MODE="$mode" nextflow run "$REPO_DIR/$harness" \
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

run_expected_failure alignment_nonzero fail \
    tests/modules/alignment/test_starsolo_fail_closed.nf \
    'STAR alignment failed'
run_expected_failure alignment_empty success \
    tests/modules/alignment/test_starsolo_fail_closed.nf \
    'STAR produced no alignment BAM'
run_expected_failure alignment_matrix_missing matrix_missing \
    tests/modules/alignment/test_starsolo_fail_closed.nf \
    'STAR produced no filtered cell matrix'
run_expected_failure alignment_bai_missing bai_missing \
    tests/modules/alignment/test_starsolo_fail_closed.nf \
    'samtools produced no BAM index'
run_expected_failure reference_nonzero fail \
    tests/modules/reference/test_prepare_reference_fail_closed.nf \
    'STAR index generation failed'
run_expected_failure reference_empty success \
    tests/modules/reference/test_prepare_reference_fail_closed.nf \
    'STAR index is incomplete'
