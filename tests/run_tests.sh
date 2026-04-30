#!/usr/bin/env bash
# tests/run_tests.sh
#
# Local CI runner for scPolASeq.
#
# Runs all test scenarios (stub and optionally real) and reports per-scenario
# pass/fail with process-count validation. Meant to be run on deepthought or
# any machine with Nextflow in PATH.
#
# Usage:
#   ./tests/run_tests.sh            # stub-run only (fast, ~1 min total)
#   ./tests/run_tests.sh --real     # also run a real conda integration test
#   ./tests/run_tests.sh --scenario full  # single scenario
#
# Exit code: 0 if all enabled scenarios pass, 1 otherwise.

set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_DIR"

# ── Colours ───────────────────────────────────────────────────────────────────
RED='\033[0;31m'; GREEN='\033[0;32m'; YELLOW='\033[1;33m'; NC='\033[0m'
pass() { echo -e "${GREEN}[PASS]${NC} $*"; }
fail() { echo -e "${RED}[FAIL]${NC} $*"; }
info() { echo -e "${YELLOW}[INFO]${NC} $*"; }

# ── Parse args ────────────────────────────────────────────────────────────────
RUN_REAL=false
TARGET_SCENARIO=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --real)        RUN_REAL=true; shift ;;
        --scenario)    TARGET_SCENARIO="$2"; shift 2 ;;
        *) echo "Unknown option: $1"; exit 1 ;;
    esac
done

# ── Scenario definitions ──────────────────────────────────────────────────────
# Format: "name|params_file|expected_stub_processes"
# expected_stub_processes is the number of "Succeeded : N" lines in the log.
declare -a SCENARIOS=(
    "full|tests/config/params_full.json|22"
    "toggles_off|tests/config/params_toggles_off.json|18"
    "from_bam|tests/config/params_from_bam.json|22"
    "feedback|tests/config/params_feedback.json|22"
    "polya_refs|tests/config/params_polya_refs.json|23"
)
# NOTE: Nextflow engine flags use single dash (-stub-run, -resume, -profile).
# Double dash (--flag) is reserved for pipeline params (params.*).

# ── Module-level tests ────────────────────────────────────────────────────────
declare -a MODULE_TESTS=(
    "pas_reference_build|tests/modules/apa/test_pas_reference_build.nf"
    "prepare_polya_references|tests/modules/reference/test_prepare_polya_references.nf"
)

# ── Helper: run one scenario ──────────────────────────────────────────────────
run_scenario() {
    local name="$1"
    local params_file="$2"
    local expected_procs="$3"
    local stub_flag="${4:--stub-run}"    # -stub-run or ""
    local work_dir="/scratch/.nf_work/test/${name}_$(date +%s)"

    local mode_label="stub"
    [[ "$stub_flag" == "" ]] && mode_label="real"

    info "Running scenario '${name}' (${mode_label})..."

    local log_file
    log_file="$(mktemp)"

    if nextflow run . \
        -profile test \
        -params-file "$params_file" \
        -work-dir "$work_dir" \
        $stub_flag \
        -ansi-log false \
        > "$log_file" 2>&1; then

        local actual_procs
        actual_procs=$(grep -c "Submitted process" "$log_file" || echo "0")

        if [[ "$mode_label" == "stub" && "$actual_procs" -ne "$expected_procs" ]]; then
            fail "Scenario '${name}' (${mode_label}): expected ${expected_procs} processes, got ${actual_procs}"
            cat "$log_file"
            rm -f "$log_file"
            return 1
        fi

        pass "Scenario '${name}' (${mode_label}): ${actual_procs} processes"
        rm -f "$log_file"
        rm -rf "$work_dir"
        return 0
    else
        fail "Scenario '${name}' (${mode_label}): Nextflow exited with error"
        grep -E "ERROR|error|WARN.*Failed" "$log_file" | head -20
        rm -f "$log_file"
        return 1
    fi
}

# ── Helper: run one module test ───────────────────────────────────────────────
run_module_test() {
    local name="$1"
    local nf_file="$2"
    local work_dir="/scratch/.nf_work/test/module_${name}_$(date +%s)"

    info "Running module test '${name}' (stub)..."

    local log_file
    log_file="$(mktemp)"

    if nextflow run "$nf_file" \
        -stub-run \
        -work-dir "$work_dir" \
        -ansi-log false \
        > "$log_file" 2>&1; then

        pass "Module test '${name}'"
        rm -f "$log_file"
        rm -rf "$work_dir"
        return 0
    else
        fail "Module test '${name}'"
        cat "$log_file"
        rm -f "$log_file"
        return 1
    fi
}

# ── Main ─────────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════════"
echo "  scPolASeq test suite"
echo "  Branch : $(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo 'unknown')"
echo "  Commit : $(git rev-parse --short HEAD 2>/dev/null || echo 'unknown')"
echo "══════════════════════════════════════════════════"
echo ""

FAILED=0

# ── 1. Workflow stub-run scenarios ────────────────────────────────────────────
for scenario_def in "${SCENARIOS[@]}"; do
    IFS='|' read -r name params_file expected_procs <<< "$scenario_def"

    # Filter to target scenario if --scenario was specified
    if [[ -n "$TARGET_SCENARIO" && "$name" != "$TARGET_SCENARIO" ]]; then
        continue
    fi

    run_scenario "$name" "$params_file" "$expected_procs" "-stub-run" || FAILED=$((FAILED + 1))
done

# ── 2. Module-level tests (stub) ─────────────────────────────────────────────
echo ""
info "Module-level tests"
for module_def in "${MODULE_TESTS[@]}"; do
    IFS='|' read -r name nf_file <<< "$module_def"
    run_module_test "$name" "$nf_file" || FAILED=$((FAILED + 1))
done

# ── 3. Real integration test (opt-in) ────────────────────────────────────────
if [[ "$RUN_REAL" == "true" ]]; then
    echo ""
    info "Real integration test (conda, full scenario)"
    run_scenario "full_real" "tests/config/params_full.json" "0" "" || FAILED=$((FAILED + 1))
fi

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "══════════════════════════════════════════════════"
if [[ "$FAILED" -eq 0 ]]; then
    pass "All tests passed"
else
    fail "${FAILED} test(s) failed"
fi
echo "══════════════════════════════════════════════════"
echo ""

exit "$FAILED"
