from typing import Dict, List, Tuple
import anndata
import numpy as np
import pandas as pd
import torch
from .preprocess import convert_to_cpm, load_markers_from_excel, prepare_marker_tensor
from .model import build_model, load_model_weights

IDX2CELLTYPE = {0: "Astrocyte", 1: "Glioblast", 2: "Radial-Glial"}

def run_astro_rf_net(
    input_anndata: anndata.AnnData,
    model_param_path: str,
    marker_gene_xlsx: str,
    n_markers: int = 50,
    case_insensitive: bool = False,
    fill_missing: float = 0.0,
    out_obs_col: str = "AstroRF-Net",
    device: str = None
) -> Tuple[anndata.AnnData, Dict[str, List[str]]]:

    markers = load_markers_from_excel(marker_gene_xlsx, top_k=n_markers)
    adata_cpm = convert_to_cpm(input_anndata)
    X_tensor, found_genes, missing_genes = prepare_marker_tensor(
        adata_cpm, markers, case_insensitive=case_insensitive, fill_missing=fill_missing
    )

    model = build_model(device=device or ("cuda" if torch.cuda.is_available() else "cpu"),
                        input_size=n_markers)
    model = load_model_weights(model, model_param_path)

    with torch.no_grad():
        y_pred = model.predict(X_tensor.numpy())  # skorch 接受 numpy

    celltype_predictions = [IDX2CELLTYPE.get(int(i), f"Unknown_{i}") for i in y_pred]
    input_anndata.obs[out_obs_col] = pd.Categorical(celltype_predictions)

    diag = {"found_genes": found_genes, "missing_genes": missing_genes}
    return input_anndata, diag
