# AstroRF-Net

AstroRF-Net is a neural network classifier that predicts **Astrocyte / Glioblast / Radial-Glial** labels from an input **AnnData** object.  
Pipeline: **log2(CPM+1) normalization → select top-50 RF markers → NN inference → write predictions to `adata.obs["AstroRF-Net"]`.**

## Install
```bash
pip install -r requirements.txt
pip install -e .
