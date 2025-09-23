# AstroRF-Net

AstroRF-Net is a neural network classifier that predicts **Astrocyte / Glioblast / Radial-Glial** labels from an input **AnnData** object.  
Pipeline: **log2(CPM+1) normalization → select top-50 RF markers → NN inference → write predictions to `adata.obs["AstroRF-Net"]`.**

## Install
```bash
pip install -r requirements.txt
pip install -i https://test.pypi.org/simple/ astrorfnet
```

## Preparation
- **[Example_adata](https://github.com/meiran5759-blip/AstroRF-Net/releases/download/v0.1.0/Gene_importance_3Celltype_233gene_mean.xlsx)**
- **[AstroRF-Net parameters](https://github.com/meiran5759-blip/AstroRF-Net/releases/download/v0.1.0/Gene_importance_3Celltype_233gene_mean.xlsx)**
- **[Gene importance](https://github.com/meiran5759-blip/AstroRF-Net/releases/download/v0.1.0/Gene_importance_3Celltype_233gene_mean.xlsx)**


## Usage
```python
import scanpy as sc
import astrorfnet
from astrorfnet import run_astro_rf_net
example_adata = sc.read_h5ad('/path/to/example_adata.h5ad')
example_adata_out, info = run_astro_rf_net(
    input_anndata=example_adata,
    model_param_path="/path/to/net.params",
    marker_gene_xlsx="/path/to/Gene_importance_3Celltype_233gene_mean.xlsx",
    n_markers=50,
    out_obs_col="AstroRF-Net"
)
```
