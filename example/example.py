import anndata
from astrorfnet import run_astro_rf_net

# 加载数据
H5AD = "data/example_adata.h5ad"
MARKERS = "data/Gene_importance_3Celltype_233gene_mean.xlsx"
WEIGHTS = "data/net.params"
OUT = "data/demo_glial_pred.h5ad"

adata = anndata.read_h5ad(H5AD)
adata_out, info = run_astro_rf_net(
    adata,
    model_param_path=WEIGHTS,
    marker_gene_xlsx=MARKERS,
    n_markers=50,
    case_insensitive=False,
    fill_missing=0.0,
    out_obs_col="AstroRF-Net"
)
print("Prediction counts:")
print(adata_out.obs["AstroRF-Net"].value_counts())
