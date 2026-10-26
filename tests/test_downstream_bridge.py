"""Offline bridge contract: synthetic files only, no notebook or network calls."""

import importlib.util
from pathlib import Path
import tempfile
import unittest


BRIDGE_PATH = Path(__file__).resolve().parents[1] / "bin" / "pbmc1k_downstream_bridge.py"
SPEC = importlib.util.spec_from_file_location("pbmc1k_downstream_bridge", BRIDGE_PATH)
bridge = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bridge)


class DownstreamBridgeTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix="downstream-bridge-")
        self.addCleanup(self.tmp.cleanup)
        root = Path(self.tmp.name)
        self.upstream = root / "upstream"
        self.package = root / "package"
        self.output = root / "output"
        for relative in (
            "labels/pbmc_1k_v3.clustered.h5ad",
            "alignment/pbmc_1k_v3.Solo.out/Gene/filtered/matrix.mtx",
            "apa/site_catalog.tsv", "apa_calls/apa_events.tsv",
            "apa_calls/pdui_usage_matrix.tsv",
            "sierra_quant/cluster/cluster_1.sierra_quant.tsv",
        ):
            path = self.upstream / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(b"synthetic fixture\n")
        for relative in ("config.py", "templates.py"):
            path = self.package / "scpolaseq_downstream" / relative
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text("# synthetic fixture\n", encoding="utf-8")

    def test_requires_each_nonempty_upstream_input(self):
        inputs = bridge.require_inputs(self.upstream, self.package)
        self.assertEqual(len(inputs), 7)
        matrix = self.upstream / "alignment/pbmc_1k_v3.Solo.out/Gene/filtered/matrix.mtx"
        matrix.write_bytes(b"")
        with self.assertRaisesRegex(ValueError, "matrix"):
            bridge.require_inputs(self.upstream, self.package)

    def test_requires_sierra_and_package_source(self):
        sierra = self.upstream / "sierra_quant/cluster/cluster_1.sierra_quant.tsv"
        sierra.unlink()
        with self.assertRaisesRegex(ValueError, "sierra_quant/cluster"):
            bridge.require_inputs(self.upstream, self.package)
        sierra.write_bytes(b"synthetic fixture\n")
        (self.package / "scpolaseq_downstream/config.py").unlink()
        with self.assertRaisesRegex(ValueError, "package_config"):
            bridge.require_inputs(self.upstream, self.package)

    def test_output_claim_never_reuses_existing_directory(self):
        bridge.claim_output(self.output)
        with self.assertRaises(FileExistsError):
            bridge.claim_output(self.output)
        (self.output / "sentinel").write_text("preserve", encoding="utf-8")
        with self.assertRaises(FileExistsError):
            bridge.claim_output(self.output)
        self.assertEqual((self.output / "sentinel").read_text(), "preserve")

    def test_source_fingerprint_changes_with_package_code(self):
        original = bridge.source_fingerprint(self.package)
        self.assertEqual(len(original), 64)
        (self.package / "scpolaseq_downstream/templates.py").write_text(
            "# changed fixture\n", encoding="utf-8")
        self.assertNotEqual(original, bridge.source_fingerprint(self.package))


if __name__ == "__main__":
    unittest.main()
