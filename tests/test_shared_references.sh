#!/usr/bin/env bash
# Verify shared reference_root resolution, fail-closed reference checks,
# STAR 2.7.4a <-> 2.7.11b index compatibility, and manifest hashing (#18).
set -euo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/scpolaseq-ref-test.XXXXXX")"
trap 'rm -rf "$TMP_ROOT"' EXIT

FIXTURE="$TMP_ROOT/ref_ssd"
mkdir -p "$FIXTURE/genomes/GRCh38" "$FIXTURE/indexes/GRCh38/star_index" "$FIXTURE/databases/polya_atlas"

cp "$REPO_DIR/tests/data/reference.fa" "$FIXTURE/genomes/GRCh38/genome.fa"
cp "$REPO_DIR/tests/data/reference.gtf" "$FIXTURE/genomes/GRCh38/genes.gtf"
cp "$REPO_DIR/tests/data/known_polya.tsv" "$FIXTURE/databases/polya_atlas/polyasite2_GRCh38.tsv"

# Derive exact contig lengths (and include CRLF to match staged GRCh38 chrom.sizes)
python3 - "$FIXTURE/genomes/GRCh38/genome.fa" "$FIXTURE/genomes/GRCh38/chrom.sizes" \
          "$FIXTURE/indexes/GRCh38/star_index/chrNameLength.txt" <<'EOF'
import sys
from pathlib import Path
sys.path.insert(0, str(Path("bin").resolve()))
from prepare_reference_bundle import read_fasta_index

entries = read_fasta_index(Path(sys.argv[1]))
with open(sys.argv[2], "wb") as handle:
    for contig, length, *_ in reversed(entries):
        handle.write(f"{contig}\t{length}\r\n".encode("utf-8"))
with open(sys.argv[3], "w", encoding="utf-8") as handle:
    for contig, length, *_ in entries:
        handle.write(f"{contig}\t{length}\n")
EOF

cat > "$FIXTURE/indexes/GRCh38/star_index/genomeParameters.txt" <<'EOF'
### STAR --runMode genomeGenerate
versionGenome	2.7.4a
genomeType	Full
sjdbOverhang	100
genomeFileSizes	152 1024
EOF
printf 'dummy_genome\n' > "$FIXTURE/indexes/GRCh38/star_index/Genome"
printf 'dummy_sa\n' > "$FIXTURE/indexes/GRCh38/star_index/SA"
printf 'dummy_saindex\n' > "$FIXTURE/indexes/GRCh38/star_index/SAindex"

run_bundle() {
    local out_dir="$1"
    shift
    mkdir -p "$out_dir"
    python3 "$REPO_DIR/bin/prepare_reference_bundle.py" "$@" \
        --out-fasta "$out_dir/ref.fa" \
        --out-gtf "$out_dir/ref.gtf" \
        --out-fai "$out_dir/ref.fa.fai" \
        --out-dict "$out_dir/ref.dict" \
        --out-sizes "$out_dir/ref.chrom.sizes" \
        --out-known-polya "$out_dir/polya.tsv" \
        --out-blacklist "$out_dir/blacklist.bed" \
        --out-manifest "$out_dir/reference_manifest.json"
}

# 1. Valid prebuilt STAR 2.7.4a index + STAR 2.7.11b + CRLF/reordered chrom.sizes
run_bundle "$TMP_ROOT/ok" \
    --genome-fasta "$FIXTURE/genomes/GRCh38/genome.fa" \
    --gtf "$FIXTURE/genomes/GRCh38/genes.gtf" \
    --known-polya "$FIXTURE/databases/polya_atlas/polyasite2_GRCh38.tsv" \
    --priming-blacklist "$REPO_DIR/assets/NO_FILE_PRIMING_BLACKLIST" \
    --chrom-sizes "$FIXTURE/genomes/GRCh38/chrom.sizes" \
    --prebuilt-star-index "$FIXTURE/indexes/GRCh38/star_index" \
    --reference-root "$FIXTURE" \
    --reference-build GRCh38 \
    --star-binary-version 2.7.11b

python3 "$REPO_DIR/bin/prepare_reference_bundle.py" \
    --finalize-star-index "$FIXTURE/indexes/GRCh38/star_index" \
    --star-index-mode prebuilt \
    --star-index-source "indexes/GRCh38/star_index" \
    --star-binary-version 2.7.11b \
    --out-manifest "$TMP_ROOT/ok/reference_manifest.json"

python3 - "$TMP_ROOT/ok/reference_manifest.json" "$FIXTURE" <<'EOF'
import json, sys
m = json.loads(open(sys.argv[1]).read())
assert m["reference_root"] == sys.argv[2]
assert m["reference_build"] == "GRCh38"
assert m["star_binary_version"] == "2.7.11b"
assert len(m["inputs"]["genome_fasta"]["sha256"]) == 64
assert len(m["inputs"]["gtf"]["sha256"]) == 64
assert len(m["inputs"]["chrom_sizes"]["sha256"]) == 64
assert m["inputs"]["known_polya"]["present"] is True
assert len(m["inputs"]["known_polya"]["sha256"]) == 64
idx = m["star_index"]
assert idx["compatible"] is True
assert idx["mode"] == "prebuilt"
assert idx["version_genome"] == "2.7.4a"
assert len(idx["parameters_sha256"]) == 64
assert len(idx["index_sha256"]) == 64
EOF
echo "PASS: valid STAR 2.7.4a index + chrom.sizes + SHA-256 manifest verified"

# 2. Fail-closed checks: empty FASTA, mismatched chrom.sizes, incompatible versionGenome, mismatched chrNameLength
expect_fail() {
    local name="$1" expected="$2"
    shift 2
    local log="$TMP_ROOT/fail_${name}.log"
    if "$@" >"$log" 2>&1; then
        echo "FAIL: ${name} unexpectedly succeeded" >&2
        return 1
    fi
    if ! grep -F "$expected" "$log" >/dev/null; then
        echo "FAIL: ${name} did not report '${expected}':" >&2
        cat "$log" >&2
        return 1
    fi
    echo "PASS: ${name} failed closed"
}

: > "$TMP_ROOT/empty.fa"
expect_fail empty_fasta "missing or empty" \
    run_bundle "$TMP_ROOT/f1" \
    --genome-fasta "$TMP_ROOT/empty.fa" \
    --gtf "$FIXTURE/genomes/GRCh38/genes.gtf" \
    --known-polya "$REPO_DIR/assets/NO_FILE" \
    --priming-blacklist "$REPO_DIR/assets/NO_FILE_PRIMING_BLACKLIST"

printf "chr_bad\t12345\n" > "$TMP_ROOT/bad.chrom.sizes"
expect_fail bad_chrom_sizes "do not match genome_fasta" \
    run_bundle "$TMP_ROOT/f2" \
    --genome-fasta "$FIXTURE/genomes/GRCh38/genome.fa" \
    --gtf "$FIXTURE/genomes/GRCh38/genes.gtf" \
    --known-polya "$REPO_DIR/assets/NO_FILE" \
    --priming-blacklist "$REPO_DIR/assets/NO_FILE_PRIMING_BLACKLIST" \
    --chrom-sizes "$TMP_ROOT/bad.chrom.sizes"

mkdir -p "$TMP_ROOT/bad_idx_ver"
cp -r "$FIXTURE/indexes/GRCh38/star_index/." "$TMP_ROOT/bad_idx_ver/"
sed -i 's/2.7.4a/2.7.1a/' "$TMP_ROOT/bad_idx_ver/genomeParameters.txt"
expect_fail bad_star_version "Incompatible STAR index versionGenome='2.7.1a'" \
    run_bundle "$TMP_ROOT/f3" \
    --genome-fasta "$FIXTURE/genomes/GRCh38/genome.fa" \
    --gtf "$FIXTURE/genomes/GRCh38/genes.gtf" \
    --known-polya "$REPO_DIR/assets/NO_FILE" \
    --priming-blacklist "$REPO_DIR/assets/NO_FILE_PRIMING_BLACKLIST" \
    --prebuilt-star-index "$TMP_ROOT/bad_idx_ver"

mkdir -p "$TMP_ROOT/bad_idx_chr"
cp -r "$FIXTURE/indexes/GRCh38/star_index/." "$TMP_ROOT/bad_idx_chr/"
printf "chr_other\t999\n" > "$TMP_ROOT/bad_idx_chr/chrNameLength.txt"
expect_fail bad_star_contigs "STAR index contigs/lengths" \
    run_bundle "$TMP_ROOT/f4" \
    --genome-fasta "$FIXTURE/genomes/GRCh38/genome.fa" \
    --gtf "$FIXTURE/genomes/GRCh38/genes.gtf" \
    --known-polya "$REPO_DIR/assets/NO_FILE" \
    --priming-blacklist "$REPO_DIR/assets/NO_FILE_PRIMING_BLACKLIST" \
    --prebuilt-star-index "$TMP_ROOT/bad_idx_chr"

# 3. End-to-end Nextflow resolution with fixture reference_root (both generic & deepthought profile override)
cat > "$TMP_ROOT/ci_resources.config" <<'EOF'
apptainer.enabled = false
conda.enabled = false
params {
    apptainer_cache_dir = null
}
process {
    cpus = 1
    memory = '2 GB'
    time = '1 h'
}
EOF

nextflow run "$REPO_DIR" \
    -profile test \
    -c "$TMP_ROOT/ci_resources.config" \
    --input "$REPO_DIR/tests/data/samplesheet.csv" \
    --reference_root "$FIXTURE" \
    --reference_build GRCh38 \
    --genome_fasta null \
    --gtf null \
    --known_polya databases/polya_atlas/polyasite2_GRCh38.tsv \
    --outdir "$TMP_ROOT/nf_generic_out" \
    -work-dir "$TMP_ROOT/nf_generic_work" \
    -stub-run -ansi-log false >/dev/null
echo "PASS: Nextflow generic --reference_root + --reference_build stub-run"

nextflow run "$REPO_DIR" \
    -profile test_pbmc1k_deepthought \
    -c "$TMP_ROOT/ci_resources.config" \
    --input "$REPO_DIR/tests/data/samplesheet.csv" \
    --reference_root "$FIXTURE" \
    --outdir "$TMP_ROOT/nf_profile_out" \
    -work-dir "$TMP_ROOT/nf_profile_work" \
    -stub-run -ansi-log false >/dev/null
echo "PASS: Nextflow test_pbmc1k_deepthought overridden with --reference_root fixture"

# 4. Real (non-stub) Nextflow PREPARE_REFERENCE_BUNDLE execution against fixture prebuilt index
nextflow run "$REPO_DIR/tests/modules/reference/test_prepare_reference_real.nf" \
    -profile test \
    -c "$TMP_ROOT/ci_resources.config" \
    --fixture "$FIXTURE" \
    --reference_root "$FIXTURE" \
    --reference_build GRCh38 \
    --outdir "$TMP_ROOT/nf_real_out" \
    -work-dir "$TMP_ROOT/nf_real_work" \
    -ansi-log false >/dev/null

python3 - "$TMP_ROOT/nf_real_out/reference/reference_manifest.json" "$FIXTURE" <<'EOF'
import json, sys
m = json.loads(open(sys.argv[1]).read())
assert m["reference_root"] == sys.argv[2]
assert m["reference_build"] == "GRCh38"
assert m["star_index"]["compatible"] is True
assert m["star_index"]["mode"] == "prebuilt"
assert m["star_index"]["version_genome"] == "2.7.4a"
assert len(m["star_index"]["index_sha256"]) == 64
assert len(m["inputs"]["genome_fasta"]["sha256"]) == 64
assert len(m["inputs"]["chrom_sizes"]["sha256"]) == 64
EOF
echo "PASS: Real Nextflow PREPARE_REFERENCE_BUNDLE validated fixture index and wrote SHA-256 manifest"
