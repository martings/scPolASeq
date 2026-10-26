#!/usr/bin/env bash
# Offline contract test: fake Apptainer, synthetic data, no image build or host data.
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TEST_DIR="$(mktemp -d "${TMPDIR:-/tmp}/scpolaseq-sif-validation.XXXXXX")"
trap 'rm -rf -- "$TEST_DIR"' EXIT
mkdir -p "$TEST_DIR/bin" "$TEST_DIR/images"
printf 'synthetic image fixture\n' > "$TEST_DIR/images/scpolaseq-python.sif"

cat > "$TEST_DIR/bin/apptainer" <<'APPTAINER'
#!/usr/bin/env bash
set -euo pipefail
[[ "$1" == exec ]] || exit 70
shift
while [[ "$1" == --* ]]; do
    case "$1" in
        --cleanenv) shift ;;
        --bind) shift 2 ;;
        *) exit 70 ;;
    esac
done
[[ "$1" == */scpolaseq-python.sif ]] || exit 70
shift
case "$1" in
    true|python) exit 0 ;;
    bedGraphToBigWig)
        case "${APPTAINER_TEST_MODE:-ok}" in
            ok) printf '\x26\xfc\x8f\x88synthetic\n' > "$4" ;;
            empty) : > "$4" ;;
            invalid) printf 'not-a-bigwig\n' > "$4" ;;
            fail) exit 23 ;;
            missing) exit 127 ;;
            *) exit 70 ;;
        esac ;;
    *) exit 70 ;;
esac
APPTAINER
chmod +x "$TEST_DIR/bin/apptainer"

run_case() {
    local mode="$1"
    PATH="$TEST_DIR/bin:$PATH" APPTAINER_TEST_MODE="$mode" \
        bash "$REPO_DIR/containers/validate.sh" "$TEST_DIR/images" --env python \
        > "$TEST_DIR/$mode.log" 2>&1
}

run_case ok || { cat "$TEST_DIR/ok.log" >&2; exit 1; }
for mode in missing empty invalid fail; do
    if run_case "$mode"; then
        echo "FAIL: Python SIF validator accepted $mode conversion" >&2
        exit 1
    fi
done
if PATH="$TEST_DIR/bin:$PATH" bash "$REPO_DIR/containers/validate.sh" \
    "$TEST_DIR/images" --env arbitrary > "$TEST_DIR/unsupported.log" 2>&1; then
    echo 'FAIL: validator accepted unsupported environment' >&2
    exit 1
fi
echo 'PASS: targeted synthetic BigWig image validation fails closed'
