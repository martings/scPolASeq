"""Check compatibility fields and fail-closed behavior of staged manifest aggregation."""

import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


HELPER = Path(__file__).resolve().parents[1] / "bin" / "scapture_manifest.py"


class ManifestTest(unittest.TestCase):
    def test_complete_manifest_and_incomplete_quant(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)

            def write(name, rows):
                path = root / name
                with path.open("w", newline="") as handle:
                    csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
                return path

            duplicate = write("duplicates.tsv", [["gene_count", "original_gene_name"], ["2", "A"]])
            pascall = write("pascall.tsv", [["metric", "value"],
                *[[f"{cat}_{kind}", "2"] for cat in ("exonic", "intronic", "3primeExtended") for kind in ("valid", "scored")],
                ["valid_peaks", "6"], ["scored_peaks", "6"]])
            reference = write("reference.tsv", [["category", "integrated_rows"],
                ["exonic", "2"], ["intronic", "0"], ["3primeExtended", "1"], ["total", "3"]])
            catalog = write("catalog.tsv", [["metric", "value"], ["site_catalog_rows", "4"]])
            quant = write("quant.tsv", [["field", "value"], ["status", "complete"],
                ["quant_matrix", "lib.KeepCell.UMIs.tsv.gz"], ["matrix_site_rows", "3"], ["matrix_cell_columns", "2"]])
            output = root / "manifest.tsv"
            args = [sys.executable, str(HELPER), "--library-id", "lib", "--bam", "lib.bam",
                "--gtf", "genes.gtf", "--fasta", "ref.fa", "--fai", "ref.fa.fai",
                "--prepared-gtf", "lib.scapture.annotation.gtf", "--gene-map", "lib.gene_name_map.tsv",
                "--duplicate-report", str(duplicate), "--pascall-metrics", str(pascall),
                "--reference-metrics", str(reference), "--catalog-metrics", str(catalog),
                "--quant-manifest", str(quant), "--species", "human", "--read-length", "90",
                "--peak-width", "400", "--threshold", "0", "--output", str(output)]
            subprocess.run(args, check=True, capture_output=True, text=True)
            with output.open(newline="") as handle:
                rows = list(csv.reader(handle, delimiter="\t"))
            self.assertEqual(rows[0], ["field", "value"])
            result = dict(rows[1:])
            self.assertEqual(result["deeppass_threshold"], "0")
            self.assertEqual(result["site_catalog_rows"], "4")
            self.assertEqual(result["merged_categories"], "2")
            self.assertEqual(result["duplicate_original_gene_names"], "1")
            self.assertEqual(result["matrix_cell_columns"], "2")
            quant = write("quant.tsv", [["field", "value"], ["status", "failed"]])
            self.assertNotEqual(subprocess.run(args, capture_output=True).returncode, 0)


if __name__ == "__main__":
    unittest.main()
