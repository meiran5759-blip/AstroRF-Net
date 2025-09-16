import argparse
import anndata
from .pipeline import run_astro_rf_net

def main():
    ap = argparse.ArgumentParser(description="AstroRF-Net inference on AnnData")
    ap.add_argument("--h5ad", required=True, help="Input AnnData (.h5ad)")
    ap.add_argument("--weights", required=True, help="Path to net.params")
    ap.add_argument("--markers", required=True, help="Excel with 'Gene' column (top 50 used)")
    ap.add_argument("--out", required=True, help="Output .h5ad with predictions written to obs")
    ap.add_argument("--col", default="AstroRF-Net", help="obs column name")
    ap.add_argument("--case-insensitive", action="store_true")
    ap.add_argument("--fill-missing", type=float, default=0.0)
    args = ap.parse_args()

    adata = anndata.read_h5ad(args.h5ad)
    adata, diag = run_astro_rf_net(
        adata,
        model_param_path=args.weights,
        marker_gene_xlsx=args.markers,
        n_markers=50,
        case_insensitive=args.case_insensitive,
        fill_missing=args.fill_missing,
        out_obs_col=args.col,
    )
    adata.write_h5ad(args.out)
    print(f"[AstroRF-Net] Done. Found {len(diag['found_genes'])}/50 markers. Missing: {len(diag['missing_genes'])}.")

