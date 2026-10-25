#!/usr/bin/env bash
# Regression tests for real BigWig output and fail-closed conversion behavior.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/scpolaseq-bigwig.XXXXXX")"
trap 'rm -rf "$TMP_ROOT"' EXIT

mkdir -p "$TMP_ROOT/data" "$TMP_ROOT/bin/no_converter" \
    "$TMP_ROOT/bin/converter_fail" "$TMP_ROOT/bin/converter_empty" \
    "$TMP_ROOT/bin/converter_ok"
printf 'synthetic bam bytes\n' > "$TMP_ROOT/data/input.bam"
printf 'synthetic bai bytes\n' > "$TMP_ROOT/data/input.bam.bai"
printf 'chr1\t1000\n' > "$TMP_ROOT/data/chrom.sizes"
: > "$TMP_ROOT/data/empty.bam"

for bin in "$TMP_ROOT"/bin/*; do
    ln -s "$(command -v awk)" "$bin/awk"
    cat > "$bin/python3" <<'PYTHON'
#!/bin/bash
if [[ "${1:-}" == */bam_to_bedgraph.py ]]; then
    shift
    while (($#)); do
        case "$1" in
            --fwd-output) fwd="$2"; shift 2 ;;
            --rev-output) rev="$2"; shift 2 ;;
            *) shift ;;
        esac
    done
    printf 'chr1\t99\t100\t1\n' > "$fwd"
    printf 'chr1\t199\t200\t2\n' > "$rev"
else
    exec /usr/bin/python3 "$@"
fi
PYTHON
    chmod +x "$bin/python3"
done

for mode in fail empty ok; do
    bin="$TMP_ROOT/bin/converter_${mode}"
    cat > "$bin/bedGraphToBigWig" <<'CONVERTER'
#!/bin/bash
case "${BIGWIG_TEST_MODE:-ok}" in
    fail) exit 23 ;;
    empty) : > "$3" ;;
    ok) printf 'synthetic nonempty BigWig bytes\n' > "$3" ;;
    *) exit 64 ;;
esac
CONVERTER
    chmod +x "$bin/bedGraphToBigWig"
    export BIGWIG_TEST_MODE="$mode"
done

run_nextflow() {
    local name="$1" bin="$2" emit="$3" bam="$4"
    local work="$TMP_ROOT/work_${name}" out="$TMP_ROOT/out_${name}"
    local log="$TMP_ROOT/${name}.log"
    local mode=ok
    [[ "$name" == converter_failure ]] && mode=fail
    [[ "$name" == empty_conversion ]] && mode=empty
    PATH="$bin:$PATH" BIGWIG_TEST_MODE="$mode" nextflow run \
        "$REPO_DIR/tests/modules/coverage/test_coverage_bigwig.nf" \
        -profile test \
        -work-dir "$work" \
        --outdir "$out" \
        --test_bam "$bam" \
        --test_bai "$TMP_ROOT/data/input.bam.bai" \
        --test_chrom_sizes "$TMP_ROOT/data/chrom.sizes" \
        --emit_bigwigs "$emit" \
        -ansi-log false > "$log" 2>&1
}

run_expected_failure() {
    local name="$1" bin="$2" emit="$3" bam="$4" expected="$5"
    local work="$TMP_ROOT/work_${name}"
    local rc=0
    if run_nextflow "$name" "$bin" "$emit" "$bam"; then
        echo "FAIL: ${name} unexpectedly succeeded" >&2
        cat "$TMP_ROOT/${name}.log" >&2
        return 1
    else
        rc=$?
    fi
    if [[ "$rc" -eq 0 ]] || ! find "$work" -type f -name '.command.err' \
        -exec grep -F "$expected" {} + >/dev/null; then
        echo "FAIL: ${name} did not record expected error '$expected' (exit ${rc})" >&2
        cat "$TMP_ROOT/${name}.log" >&2
        find "$work" -name '.command.err' -type f -print -exec cat {} \; >&2 || true
        return 1
    fi
    echo "PASS: ${name} failed closed as expected"
}

run_expected_failure missing_converter "$TMP_ROOT/bin/no_converter" true \
    "$TMP_ROOT/data/input.bam" 'bedGraphToBigWig executable unavailable'
run_expected_failure converter_failure "$TMP_ROOT/bin/converter_fail" true \
    "$TMP_ROOT/data/input.bam" 'Forward bedGraphToBigWig conversion failed'
run_expected_failure empty_conversion "$TMP_ROOT/bin/converter_empty" true \
    "$TMP_ROOT/data/input.bam" 'Forward BigWig output is missing or empty'
run_expected_failure empty_bam "$TMP_ROOT/bin/converter_ok" false \
    "$TMP_ROOT/data/empty.bam" 'Coverage input BAM is missing or empty'

if ! run_nextflow successful_conversion "$TMP_ROOT/bin/converter_ok" true \
    "$TMP_ROOT/data/input.bam"; then
    cat "$TMP_ROOT/successful_conversion.log" >&2
    exit 1
fi
test -s "$TMP_ROOT/out_successful_conversion/coverage/cluster/cluster_1.fwd.bw"
test -s "$TMP_ROOT/out_successful_conversion/coverage/cluster/cluster_1.rev.bw"
echo "PASS: successful_conversion produced two nonempty BigWig files"

if ! run_nextflow bigwigs_disabled "$TMP_ROOT/bin/no_converter" false \
    "$TMP_ROOT/data/input.bam"; then
    cat "$TMP_ROOT/bigwigs_disabled.log" >&2
    exit 1
fi
test -s "$TMP_ROOT/out_bigwigs_disabled/coverage/cluster/cluster_1.fwd.bedGraph"
if find "$TMP_ROOT/out_bigwigs_disabled" -type f -name '*.bw' | grep -q .; then
    echo 'FAIL: BigWig outputs were created while disabled' >&2
    exit 1
fi
echo "PASS: bigwigs_disabled produced bedGraphs and no placeholder BigWigs"
