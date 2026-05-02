import csv
import importlib.util
import subprocess
import sys
import tempfile
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]


def load_module(relative_path: str, module_name: str):
    module_path = REPO_ROOT / relative_path
    spec = importlib.util.spec_from_file_location(module_name, module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def write_tsv(path: Path, fieldnames, rows) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def read_tsv(path: Path):
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class DummyAdata:
    def __init__(self, obs_names):
        self.obs_names = list(obs_names)
        self.obs = {}
        self.obsm = {}

    @property
    def n_obs(self):
        return len(self.obs_names)


def test_extract_barcode_registry_preserves_starsolo_suffixes(tmp_path):
    extract = load_module("bin/extract_barcode_registry.py", "extract_barcode_registry")
    solo_dir = tmp_path / "Solo.out" / "Gene" / "filtered"
    solo_dir.mkdir(parents=True)
    (solo_dir / "barcodes.tsv").write_text("AAACCCAAGGAGAGTA-1\nTTTGTTTCCCAAATCG-2\n", encoding="utf-8")

    rows = extract.extract_from_solo(tmp_path / "Solo.out", "pbmc_1k_v3", "pbmc_1k_v3_L002")

    assert [row["barcode_raw"] for row in rows] == ["AAACCCAAGGAGAGTA-1", "TTTGTTTCCCAAATCG-2"]
    assert [row["barcode_corrected"] for row in rows] == ["AAACCCAAGGAGAGTA-1", "TTTGTTTCCCAAATCG-2"]
    assert {row["library_id"] for row in rows} == {"pbmc_1k_v3_L002"}


def test_scanpy_cluster_outputs_library_aware_annotations_and_reversible_ids(tmp_path):
    scanpy_cluster = load_module("bin/scanpy_cluster_sc.py", "scanpy_cluster_sc")

    matrix_dir = tmp_path / "matrix" / "filtered"
    matrix_dir.mkdir(parents=True)
    (matrix_dir / "barcodes.tsv").write_text("AAAC-1\nTTTG-1\n", encoding="utf-8")

    metadata_path = tmp_path / "cell_metadata.tsv"
    metadata_path.write_text("", encoding="utf-8")
    out_annotations = tmp_path / "annotations.tsv"
    out_report = tmp_path / "report.tsv"
    out_h5ad = tmp_path / "clustered.h5ad"

    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "bin/scanpy_cluster_sc.py"),
            "--matrix-dir",
            str(tmp_path / "matrix"),
            "--cell-metadata",
            str(metadata_path),
            "--sample-id",
            "pbmc_1k_v3",
            "--library-id",
            "pbmc_1k_v3_L002",
            "--enable-internal-clustering",
            "false",
            "--out-annotations",
            str(out_annotations),
            "--out-report",
            str(out_report),
            "--out-h5ad",
            str(out_h5ad),
        ],
        check=True,
        cwd=REPO_ROOT,
    )

    rows = read_tsv(out_annotations)
    assert rows[0]["library_id"] == "pbmc_1k_v3_L002"
    assert rows[0]["barcode_raw"] == "AAAC-1"
    assert rows[0]["barcode_corrected"] == "AAAC-1"
    assert rows[0]["cell_id"] == "pbmc_1k_v3:pbmc_1k_v3_L002:AAAC-1"

    dummy = DummyAdata(["AAAC-1", "TTTG-1"])
    scanpy_cluster.annotate_adata_obs(dummy, "pbmc_1k_v3", "pbmc_1k_v3_L002")
    assert dummy.obs["sample_id"] == ["pbmc_1k_v3", "pbmc_1k_v3"]
    assert dummy.obs["library_id"] == ["pbmc_1k_v3_L002", "pbmc_1k_v3_L002"]
    assert dummy.obs["barcode_raw"] == ["AAAC-1", "TTTG-1"]
    assert dummy.obs["barcode_corrected"] == ["AAAC-1", "TTTG-1"]
    assert dummy.obs["cell_id"] == [
        "pbmc_1k_v3:pbmc_1k_v3_L002:AAAC-1",
        "pbmc_1k_v3:pbmc_1k_v3_L002:TTTG-1",
    ]
    assert dummy.obs_names == dummy.obs["cell_id"]


def test_build_group_map_keeps_library_id_for_shared_barcodes(tmp_path):
    annotation_fields = [
        "sample_id",
        "library_id",
        "barcode_raw",
        "barcode_corrected",
        "cell_id",
        "cluster_id",
        "cell_type",
        "condition",
        "batch",
        "label_source",
    ]
    annotation_one = tmp_path / "L001.cell_annotations.tsv"
    annotation_two = tmp_path / "L002.cell_annotations.tsv"
    write_tsv(
        annotation_one,
        annotation_fields,
        [
            {
                "sample_id": "pbmc_1k_v3",
                "library_id": "pbmc_1k_v3_L001",
                "barcode_raw": "AAAC-1",
                "barcode_corrected": "AAAC-1",
                "cell_id": "pbmc_1k_v3:pbmc_1k_v3_L001:AAAC-1",
                "cluster_id": "cluster_1",
                "cell_type": "T",
                "condition": "",
                "batch": "",
                "label_source": "unit_test",
            }
        ],
    )
    write_tsv(
        annotation_two,
        annotation_fields,
        [
            {
                "sample_id": "pbmc_1k_v3",
                "library_id": "pbmc_1k_v3_L002",
                "barcode_raw": "AAAC-1",
                "barcode_corrected": "AAAC-1",
                "cell_id": "pbmc_1k_v3:pbmc_1k_v3_L002:AAAC-1",
                "cluster_id": "cluster_2",
                "cell_type": "B",
                "condition": "",
                "batch": "",
                "label_source": "unit_test",
            }
        ],
    )

    out_annotations = tmp_path / "cell_annotations.tsv"
    out_group_map = tmp_path / "group_map.tsv"
    out_group_summary = tmp_path / "group_summary.tsv"
    subprocess.run(
        [
            sys.executable,
            str(REPO_ROOT / "bin/build_group_map.py"),
            "--annotation-files",
            str(annotation_one),
            str(annotation_two),
            "--grouping-levels",
            "cell,cluster",
            "--out-cell-annotations",
            str(out_annotations),
            "--out-group-map",
            str(out_group_map),
            "--out-group-summary",
            str(out_group_summary),
        ],
        check=True,
        cwd=REPO_ROOT,
    )

    annotation_rows = read_tsv(out_annotations)
    assert set(annotation_rows[0].keys()) == set(annotation_fields)
    assert {row["library_id"] for row in annotation_rows} == {"pbmc_1k_v3_L001", "pbmc_1k_v3_L002"}

    group_rows = read_tsv(out_group_map)
    assert set(group_rows[0].keys()) == {"sample_id", "library_id", "barcode_corrected", "group_level", "group_id"}
    assert {
        (row["library_id"], row["group_level"], row["group_id"])
        for row in group_rows
    } == {
        ("pbmc_1k_v3_L001", "cell", "pbmc_1k_v3:pbmc_1k_v3_L001:AAAC-1"),
        ("pbmc_1k_v3_L001", "cluster", "cluster_1"),
        ("pbmc_1k_v3_L002", "cell", "pbmc_1k_v3:pbmc_1k_v3_L002:AAAC-1"),
        ("pbmc_1k_v3_L002", "cluster", "cluster_2"),
    }


def test_grouped_3prime_coverage_scopes_shared_barcodes_by_library(tmp_path):
    grouped_coverage = load_module("bin/grouped_3prime_coverage.py", "grouped_3prime_coverage")
    group_map_path = tmp_path / "group_map.tsv"
    barcode_registry_path = tmp_path / "barcode_registry.tsv"

    write_tsv(
        group_map_path,
        ["sample_id", "library_id", "barcode_corrected", "group_level", "group_id"],
        [
            {
                "sample_id": "pbmc_1k_v3",
                "library_id": "pbmc_1k_v3_L001",
                "barcode_corrected": "AAAC-1",
                "group_level": "cluster",
                "group_id": "cluster_1",
            },
            {
                "sample_id": "pbmc_1k_v3",
                "library_id": "pbmc_1k_v3_L002",
                "barcode_corrected": "AAAC-1",
                "group_level": "cluster",
                "group_id": "cluster_2",
            },
        ],
    )
    write_tsv(
        barcode_registry_path,
        ["sample_id", "library_id", "barcode_raw", "barcode_corrected", "source"],
        [
            {
                "sample_id": "pbmc_1k_v3",
                "library_id": "pbmc_1k_v3_L001",
                "barcode_raw": "AAAC-1",
                "barcode_corrected": "AAAC-1",
                "source": "starsolo",
            },
            {
                "sample_id": "pbmc_1k_v3",
                "library_id": "pbmc_1k_v3_L002",
                "barcode_raw": "AAAC-1",
                "barcode_corrected": "AAAC-1",
                "source": "starsolo",
            },
        ],
    )

    group_map = grouped_coverage.load_group_map(group_map_path, "pbmc_1k_v3", "pbmc_1k_v3_L001")
    barcode_registry = grouped_coverage.load_barcode_registry(
        barcode_registry_path,
        "pbmc_1k_v3",
        "pbmc_1k_v3_L001",
    )
    counts = grouped_coverage.fallback_counts(group_map, barcode_registry, [("chr1", 1000)])

    assert group_map == {"AAAC-1": [("cluster", "cluster_1")]}
    assert len(barcode_registry) == 1
    assert {key[1] for key in counts} == {"cluster_1"}


def test_scanpy_cluster_emits_embedding_plots_when_embeddings_exist(tmp_path):
    scanpy_cluster = load_module("bin/scanpy_cluster_sc.py", "scanpy_cluster_sc_plots")
    if importlib.util.find_spec("matplotlib") is None:
        return

    import numpy as np

    dummy = DummyAdata(["AAAC-1", "TTTG-1"])
    dummy.obsm["X_pca"] = np.array([[0.0, 1.0], [1.0, 0.0]])
    dummy.obsm["X_umap"] = np.array([[0.2, 0.3], [0.8, 0.5]])
    dummy.obs["cluster_id"] = ["cluster_1", "cluster_2"]

    plot_prefix = tmp_path / "pbmc_1k_v3_L001.clustered"
    scanpy_cluster.emit_embedding_plots(dummy, str(plot_prefix))

    assert (tmp_path / "pbmc_1k_v3_L001.clustered.pca.png").exists()
    assert (tmp_path / "pbmc_1k_v3_L001.clustered.umap.png").exists()


def test_write_group_bams_names_outputs_with_library_prefix(tmp_path):
    write_group_bams = load_module("bin/write_group_bams.py", "write_group_bams")
    assert (
        write_group_bams.grouped_bam_name("pbmc_1k_v3_L001", "cluster_1")
        == "pbmc_1k_v3_L001.cluster_1.grouped.bam"
    )
    assert write_group_bams.grouped_bam_name(None, "cluster_1") == "cluster_1.grouped.bam"


def test_alignment_subworkflow_collapses_samplesheet_rows_by_sample_id(tmp_path):
    del tmp_path
    alignment_nf = (REPO_ROOT / "subworkflows" / "local" / "alignment_or_import" / "main.nf").read_text(
        encoding="utf-8"
    )
    assert ".groupTuple()" in alignment_nf
    assert "library_id        : sampleId" in alignment_nf
    assert "source_library_ids" in alignment_nf


def test_apa_core_normalizes_group_id_before_sierra(tmp_path):
    del tmp_path
    apa_core_nf = (REPO_ROOT / "subworkflows" / "apa_core.nf").read_text(encoding="utf-8")
    assert "def prefix = meta.library_id ? \"${meta.library_id}.\" : ''" in apa_core_nf
    assert "group_id.startsWith(prefix)" in apa_core_nf


def run_all_tests() -> None:
    test_functions = [
        test_extract_barcode_registry_preserves_starsolo_suffixes,
        test_scanpy_cluster_outputs_library_aware_annotations_and_reversible_ids,
        test_build_group_map_keeps_library_id_for_shared_barcodes,
        test_grouped_3prime_coverage_scopes_shared_barcodes_by_library,
        test_scanpy_cluster_emits_embedding_plots_when_embeddings_exist,
        test_write_group_bams_names_outputs_with_library_prefix,
        test_alignment_subworkflow_collapses_samplesheet_rows_by_sample_id,
        test_apa_core_normalizes_group_id_before_sierra,
    ]
    for test_function in test_functions:
        with tempfile.TemporaryDirectory() as tmpdir:
            test_function(Path(tmpdir))
    print(f"{len(test_functions)} tests passed")


if __name__ == "__main__":
    run_all_tests()
