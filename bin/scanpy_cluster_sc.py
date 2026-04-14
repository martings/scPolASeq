#!/usr/bin/env python3
"""
scanpy_cluster_sc.py — Cell annotation and clustering for scPolASeq.

Label source priority (first available wins):
  1. External TSV metadata          (--cell-metadata)
  2. CellTypist auto-annotation     (--celltypist-model)   [requires celltypist pkg]
  3. h5ad reference label transfer  (--reference-h5ad  --reference-label-col)
  4. Scanpy internal clustering     (--enable-internal-clustering true)
     - Leiden via leidenalg+igraph  (if installed)
     - KMeans on PCA coordinates    (fallback, always available)
  5. Round-robin fallback           (last resort)
"""

import argparse
import csv
import gzip
import json
import logging
import warnings
from pathlib import Path

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
log = logging.getLogger()

CANONICAL_COLUMNS = [
    "sample_id", "barcode_raw", "barcode_corrected", "cell_id",
    "cluster_id", "cell_type", "condition", "batch", "label_source",
]

_SENTINEL = ("None", "NO_FILE", "null", "")


# ── helpers ───────────────────────────────────────────────────────────────────

def _is_sentinel(value: str) -> bool:
    return not value or value.strip() in _SENTINEL


def load_metadata(path: Path, sample_id: str):
    if not path.exists() or _is_sentinel(path.name) or path.stat().st_size == 0:
        return []
    with path.open("r", newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [r for r in reader if (r.get("sample_id") or "") == sample_id]


def _pick_barcode_files(matrix_dir: Path):
    """Prefer STARsolo filtered/ barcodes over raw/ to avoid empty-droplet inflation."""
    for subdir in ("filtered", "raw", ""):
        pattern = f"{subdir}/barcodes.tsv" if subdir else "barcodes.tsv"
        hits = list(matrix_dir.rglob(pattern)) + list(matrix_dir.rglob(pattern + ".gz"))
        if hits:
            return hits
    return []


def find_barcodes(matrix_dir: Path):
    barcodes = []
    for candidate in _pick_barcode_files(matrix_dir):
        opener = gzip.open(candidate, "rt", encoding="utf-8") if candidate.suffix == ".gz" \
            else candidate.open("r", encoding="utf-8")
        with opener as fh:
            for line in fh:
                bc = line.strip()
                if bc:
                    barcodes.append(bc)
    return barcodes


def find_filtered_matrix_dir(matrix_dir: Path) -> Path:
    """Return the 10x-format subdirectory containing matrix.mtx[.gz]."""
    for subdir in ("Gene/filtered", "GeneFull/filtered", "filtered", ""):
        candidate = matrix_dir / subdir if subdir else matrix_dir
        if (candidate / "matrix.mtx").exists() or (candidate / "matrix.mtx.gz").exists():
            return candidate
    return matrix_dir


def _load_10x_matrix(mtx_dir: Path):
    """Load a 10x-format directory into AnnData, handling both gzipped and plain files."""
    import scipy.io
    import anndata
    import gzip

    def _open(path_plain, path_gz, mode="r"):
        if path_gz.exists():
            return gzip.open(str(path_gz), mode + "t", encoding="utf-8")
        return path_plain.open(mode, encoding="utf-8")

    mtx_plain = mtx_dir / "matrix.mtx"
    mtx_gz    = mtx_dir / "matrix.mtx.gz"
    bc_plain  = mtx_dir / "barcodes.tsv"
    bc_gz     = mtx_dir / "barcodes.tsv.gz"

    # Features file: STARsolo uses features.tsv, Cell Ranger uses genes.tsv
    ft_plain = mtx_dir / "features.tsv" if (mtx_dir / "features.tsv").exists() else mtx_dir / "genes.tsv"
    ft_gz    = Path(str(ft_plain) + ".gz")

    if mtx_gz.exists():
        with gzip.open(str(mtx_gz), "rb") as fh:
            mat = scipy.io.mmread(fh).T.tocsr()
    elif mtx_plain.exists():
        mat = scipy.io.mmread(str(mtx_plain)).T.tocsr()
    else:
        raise FileNotFoundError(f"No matrix.mtx[.gz] in {mtx_dir}")

    with _open(bc_plain, bc_gz) as fh:
        barcodes = [line.strip() for line in fh if line.strip()]

    with _open(ft_plain, ft_gz) as fh:
        # features.tsv: gene_id<tab>gene_name<tab>feature_type
        var_names = []
        for line in fh:
            parts = line.strip().split("\t")
            var_names.append(parts[1] if len(parts) > 1 else parts[0])

    adata = anndata.AnnData(X=mat)
    adata.obs_names = barcodes
    adata.var_names = var_names
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", UserWarning)
        adata.var_names_make_unique()
    return adata


def _make_row(sample_id, bc, cluster_id, cell_type, label_source, condition="", batch=""):
    return {
        "sample_id": sample_id, "barcode_raw": bc, "barcode_corrected": bc,
        "cell_id": f"{sample_id}:{bc}", "cluster_id": cluster_id,
        "cell_type": cell_type, "condition": condition, "batch": batch,
        "label_source": label_source,
    }


# ── label sources ─────────────────────────────────────────────────────────────

def labels_from_celltypist(matrix_dir: Path, model_arg: str, barcodes: list, sample_id: str):
    """Run CellTypist. Returns (rows, mode) or (None, None) on any failure."""
    try:
        import celltypist
        import scanpy as sc
    except ImportError as exc:
        log.warning(f"celltypist/scanpy not importable ({exc}); skipping CellTypist tier.")
        return None, None

    mtx_dir = find_filtered_matrix_dir(matrix_dir)
    log.info(f"CellTypist: loading matrix from {mtx_dir}")
    try:
        adata = _load_10x_matrix(mtx_dir)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        model = celltypist.models.Model.load(model_arg)
        pred = celltypist.annotate(adata, model=model, majority_voting=True)
        pred_adata = pred.to_adata()
        ct_map = pred_adata.obs["majority_voting"].to_dict()
    except Exception as exc:
        log.warning(f"CellTypist prediction failed: {exc}; skipping.")
        return None, None

    rows = [_make_row(sample_id, bc, ct_map.get(bc, "unknown"), ct_map.get(bc, "unknown"),
                      "celltypist") for bc in barcodes]
    n_types = len(set(r["cell_type"] for r in rows))
    log.info(f"CellTypist: {n_types} cell types across {len(rows)} cells")
    return rows, "celltypist"


def labels_from_reference_h5ad(ref_path: Path, label_col: str, barcodes: list, sample_id: str):
    """Transfer labels from reference h5ad by barcode matching. Returns (rows, mode) or (None, None)."""
    try:
        import anndata
    except ImportError as exc:
        log.warning(f"anndata not importable ({exc}); skipping h5ad reference tier.")
        return None, None

    try:
        ref = anndata.read_h5ad(str(ref_path))
    except Exception as exc:
        log.warning(f"Failed to load reference h5ad {ref_path}: {exc}; skipping.")
        return None, None

    if label_col not in ref.obs.columns:
        available = list(ref.obs.columns)
        log.warning(f"Column '{label_col}' not in reference h5ad obs (available: {available}); skipping.")
        return None, None

    bc_to_label = {str(idx): str(ref.obs.at[idx, label_col]) for idx in ref.obs.index}
    matched = sum(1 for bc in barcodes if bc in bc_to_label)
    log.info(f"h5ad reference: {matched}/{len(barcodes)} barcodes matched (col={label_col})")
    if matched == 0:
        log.warning("h5ad reference: 0 barcodes matched; skipping this tier.")
        return None, None

    rows = [_make_row(sample_id, bc, bc_to_label.get(bc, "unmatched"),
                      bc_to_label.get(bc, "unmatched"), "reference_h5ad") for bc in barcodes]
    return rows, "reference_h5ad"


def labels_from_leiden(matrix_dir: Path, barcodes: list, sample_id: str, resolution=0.5):
    """Leiden clustering via scanpy (requires leidenalg + igraph). Returns (rows, mode, adata) or (None, None, None)."""
    try:
        import leidenalg  # noqa: F401 — just to check availability
        import scanpy as sc
    except ImportError as exc:
        log.info(f"leidenalg/scanpy not available ({exc}); will try KMeans instead.")
        return None, None, None

    mtx_dir = find_filtered_matrix_dir(matrix_dir)
    log.info(f"Leiden: loading matrix from {mtx_dir}")
    try:
        adata = _load_10x_matrix(mtx_dir)
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
        adata_hvg = adata[:, adata.var["highly_variable"]].copy()
        adata_hvg.X = adata_hvg.X.toarray()
        sc.pp.scale(adata_hvg, max_value=10)
        n_comps = min(30, adata_hvg.n_obs - 1, adata_hvg.n_vars - 1)
        sc.tl.pca(adata_hvg, n_comps=n_comps)
        sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=n_comps)
        sc.tl.leiden(adata_hvg, resolution=resolution)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.tl.umap(adata_hvg)
    except Exception as exc:
        log.warning(f"Leiden clustering failed: {exc}; will try KMeans instead.")
        return None, None, None

    cluster_map = adata_hvg.obs["leiden"].to_dict()
    log.info(f"Leiden: {len(set(cluster_map.values()))} clusters for {len(cluster_map)} cells")
    rows = []
    for bc in barcodes:
        cid = f"cluster_{int(cluster_map[bc]) + 1}" if bc in cluster_map else "unassigned"
        rows.append(_make_row(sample_id, bc, cid, "unlabeled", "leiden_clustering"))
    adata_hvg.obs["cluster_id"] = [f"cluster_{int(cluster_map.get(bc, -1)) + 1}" if bc in cluster_map else "unassigned" for bc in adata_hvg.obs_names]
    adata_hvg.obs["label_source"] = "leiden_clustering"
    return rows, "leiden_clustering", adata_hvg


def labels_from_kmeans(matrix_dir: Path, barcodes: list, sample_id: str, n_clusters: int = 8):
    """KMeans on PCA coordinates via scanpy + sklearn. Returns (rows, mode, adata) or (None, None, None)."""
    try:
        import scanpy as sc
        from sklearn.cluster import KMeans
    except ImportError as exc:
        log.warning(f"scanpy/sklearn not importable ({exc}); skipping KMeans tier.")
        return None, None, None

    mtx_dir = find_filtered_matrix_dir(matrix_dir)
    log.info(f"KMeans: loading matrix from {mtx_dir}")
    try:
        adata = _load_10x_matrix(mtx_dir)
        sc.pp.filter_cells(adata, min_genes=200)
        sc.pp.filter_genes(adata, min_cells=3)
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="seurat")
        adata_hvg = adata[:, adata.var["highly_variable"]].copy()
        adata_hvg.X = adata_hvg.X.toarray()  # 2000 HVGs × n_cells — safe to densify
        sc.pp.scale(adata_hvg, max_value=10)
        n_comps = min(30, adata_hvg.n_obs - 1, adata_hvg.n_vars - 1)
        sc.tl.pca(adata_hvg, n_comps=n_comps)
    except Exception as exc:
        log.warning(f"Matrix loading/PCA failed: {exc}; skipping KMeans tier.")
        return None, None, None

    k = min(n_clusters, adata_hvg.n_obs)
    km = KMeans(n_clusters=k, random_state=42, n_init=10)
    labels = km.fit_predict(adata_hvg.obsm["X_pca"])
    cluster_map = dict(zip(adata_hvg.obs_names.tolist(), (int(l) for l in labels)))
    log.info(f"KMeans: {k} clusters for {len(cluster_map)} cells")
    rows = []
    for bc in barcodes:
        cid = f"cluster_{cluster_map[bc] + 1}" if bc in cluster_map else "unassigned"
        rows.append(_make_row(sample_id, bc, cid, "unlabeled", "kmeans_clustering"))

    # Annotate adata_hvg with cluster labels and compute UMAP for h5ad output
    adata_hvg.obs["cluster_id"] = [f"cluster_{cluster_map[bc] + 1}" if bc in cluster_map else "unassigned" for bc in adata_hvg.obs_names]
    adata_hvg.obs["label_source"] = "kmeans_clustering"
    try:
        sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=adata_hvg.obsm["X_pca"].shape[1])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc.tl.umap(adata_hvg)
    except Exception as exc:
        log.warning(f"UMAP computation failed ({exc}); h5ad will have PCA only.")
    return rows, "kmeans_clustering", adata_hvg


def fallback_round_robin(barcodes: list, sample_id: str):
    rows = [_make_row(sample_id, bc, f"cluster_{((i) % 3) + 1}", "unlabeled",
                      "fallback_internal_clustering")
            for i, bc in enumerate(barcodes)]
    return rows, "fallback_internal_clustering"


# ── main ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix-dir",               required=True)
    parser.add_argument("--cell-metadata",             required=True)
    parser.add_argument("--sample-id",                 required=True)
    parser.add_argument("--library-id",                required=True)
    parser.add_argument("--enable-internal-clustering", required=True)
    parser.add_argument("--celltypist-model",          default=None,
                        help="Path to .pkl or built-in model name e.g. 'Immune_All_Low.pkl'")
    parser.add_argument("--reference-h5ad",            default=None,
                        help="Reference h5ad for barcode-matched label transfer")
    parser.add_argument("--reference-label-col",       default="cell_type",
                        help="obs column in reference h5ad to use for labels")
    parser.add_argument("--out-annotations",           required=True)
    parser.add_argument("--out-report",                required=True)
    parser.add_argument("--out-h5ad",                  required=True)
    args = parser.parse_args()

    sample_id = args.sample_id
    matrix_dir = Path(args.matrix_dir)
    enable_clustering = args.enable_internal_clustering.lower() == "true"
    rows = None
    clustering_mode = None
    adata_out = None

    # ── Tier 1: External metadata TSV ─────────────────────────────────────────
    meta_rows = load_metadata(Path(args.cell_metadata), sample_id)
    if meta_rows:
        rows = [{
            "sample_id":         r.get("sample_id", ""),
            "barcode_raw":       r.get("barcode_raw") or r.get("barcode_corrected", ""),
            "barcode_corrected": r.get("barcode_corrected") or r.get("barcode_raw", ""),
            "cell_id":           r.get("cell_id") or f"{sample_id}:{r.get('barcode_corrected') or r.get('barcode_raw', '')}",
            "cluster_id":        r.get("cluster_id", ""),
            "cell_type":         r.get("cell_type", ""),
            "condition":         r.get("condition", ""),
            "batch":             r.get("batch", ""),
            "label_source":      r.get("label_source", "external_metadata"),
        } for r in meta_rows]
        clustering_mode = "external_metadata"
        log.info(f"Tier 1 (external metadata): {len(rows)} cells")

    # Load barcodes (needed for all internal tiers)
    barcodes = find_barcodes(matrix_dir)
    log.info(f"Barcodes from matrix: {len(barcodes)}")

    # ── Tier 2: CellTypist ────────────────────────────────────────────────────
    ct_model = args.celltypist_model
    if rows is None and ct_model and not _is_sentinel(ct_model):
        rows, clustering_mode = labels_from_celltypist(matrix_dir, ct_model, barcodes, sample_id)
        if rows is not None:
            log.info(f"Tier 2 (CellTypist): {len(rows)} cells")

    # ── Tier 3: h5ad reference label transfer ─────────────────────────────────
    ref_h5ad = args.reference_h5ad
    if rows is None and ref_h5ad and not _is_sentinel(ref_h5ad):
        ref_path = Path(ref_h5ad)
        if ref_path.exists() and ref_path.stat().st_size > 0:
            rows, clustering_mode = labels_from_reference_h5ad(
                ref_path, args.reference_label_col, barcodes, sample_id)
            if rows is not None:
                log.info(f"Tier 3 (h5ad reference): {len(rows)} cells")

    # ── Tier 4: Internal Scanpy clustering ────────────────────────────────────
    if rows is None and enable_clustering and barcodes:
        # Try Leiden first, then KMeans
        rows, clustering_mode, adata_out = labels_from_leiden(matrix_dir, barcodes, sample_id)
        if rows is None:
            rows, clustering_mode, adata_out = labels_from_kmeans(matrix_dir, barcodes, sample_id)
        if rows is not None:
            log.info(f"Tier 4 (internal clustering/{clustering_mode}): {len(rows)} cells")

    # ── Tier 5: Round-robin fallback ──────────────────────────────────────────
    if rows is None:
        if not barcodes:
            rows, clustering_mode = [], "empty_annotations"
        else:
            rows, clustering_mode = fallback_round_robin(barcodes, sample_id)
            log.info(f"Tier 5 (round-robin fallback): {len(rows)} cells")

    # ── Write outputs ─────────────────────────────────────────────────────────
    with Path(args.out_annotations).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=CANONICAL_COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    with Path(args.out_report).open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["sample_id", "library_id", "n_cells", "clustering_mode"])
        writer.writerow([sample_id, args.library_id, len(rows), clustering_mode])

    if adata_out is not None:
        adata_out.obs["sample_id"]  = sample_id
        adata_out.obs["library_id"] = args.library_id
        try:
            adata_out.write_h5ad(args.out_h5ad)
            log.info(f"Saved h5ad: {args.out_h5ad} ({adata_out.n_obs} cells x {adata_out.n_vars} genes)")
        except Exception as exc:
            log.warning(f"h5ad write failed ({exc}); writing JSON stub instead.")
            Path(args.out_h5ad).write_text(
                json.dumps({"sample_id": sample_id, "library_id": args.library_id,
                            "n_cells": len(rows), "clustering_mode": clustering_mode}, indent=2),
                encoding="utf-8")
    else:
        Path(args.out_h5ad).write_text(
            json.dumps({"sample_id": sample_id, "library_id": args.library_id,
                        "n_cells": len(rows), "clustering_mode": clustering_mode}, indent=2),
            encoding="utf-8",
        )
    log.info(f"Done: {len(rows)} cells, mode={clustering_mode}")


if __name__ == "__main__":
    main()
