from typing import List, Tuple, Dict
import numpy as np
import pandas as pd
import anndata
from scipy.sparse import issparse
import torch

def convert_to_cpm(adata: anndata.AnnData) -> anndata.AnnData:
    counts_matrix = adata.X
    if issparse(counts_matrix):
        counts_matrix = counts_matrix.toarray()
    counts_matrix = counts_matrix.astype(float, copy=False)
    total_counts_per_cell = counts_matrix.sum(axis=1)
    total_counts_per_cell = np.maximum(total_counts_per_cell, 1.0)
    cpm_matrix = (counts_matrix / total_counts_per_cell[:, None]) * 1e6
    log2_cpm_matrix = np.log2(cpm_matrix + 1.0).astype(np.float32)
    out = anndata.AnnData(X=log2_cpm_matrix, obs=adata.obs.copy(), var=adata.var.copy())
    out.var_names = adata.var_names.copy()
    out.obs_names = adata.obs_names.copy()
    return out

def load_markers_from_excel(xlsx_path: str, top_k: int = 50) -> List[str]:
    df = pd.read_excel(xlsx_path, header=0, index_col=None)
    if "Gene" not in df.columns:
        raise ValueError(f"Expected column 'Gene' in {xlsx_path}")
    genes = df["Gene"].astype(str).values[:top_k].tolist()
    if len(genes) < top_k:
        raise ValueError(f"Only {len(genes)} genes found, need {top_k}.")
    return genes

def prepare_marker_tensor(
    adata: anndata.AnnData,
    marker_genes: List[str],
    case_insensitive: bool = False,
    fill_missing: float = 0.0
) -> Tuple[torch.Tensor, List[str], List[str]]:
    var_names = adata.var_names.astype(str)
    if case_insensitive:
        var_lower = pd.Index([g.lower() for g in var_names])
        idx_map = {g: i for i, g in enumerate(var_lower)}
        idxs, found, missing = [], [], []
        for m in marker_genes:
            key = m.lower()
            if key in idx_map:
                idxs.append(idx_map[key]); found.append(m)
            else:
                idxs.append(None); missing.append(m)
    else:
        idx_map = {g: i for i, g in enumerate(var_names)}
        idxs, found, missing = [], [], []
        for m in marker_genes:
            if m in idx_map:
                idxs.append(idx_map[m]); found.append(m)
            else:
                idxs.append(None); missing.append(m)

    X = adata.X
    if issparse(X):
        X = X.toarray()
    X = X.astype(np.float32, copy=False)

    n_cells = X.shape[0]
    feat = np.empty((n_cells, len(marker_genes)), dtype=np.float32)
    for j, idx in enumerate(idxs):
        if idx is None:
            feat[:, j] = fill_missing
        else:
            feat[:, j] = X[:, idx]
    return torch.from_numpy(feat), found, missing
