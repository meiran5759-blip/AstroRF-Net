import anndata
from astrorfnet import AstroRFNet

# 加载数据
adata = anndata.read_h5ad("data/example_adata.h5ad")
marker_genes = pd.read_excel("data/marker_genes.xlsx")['Gene'].values

# 初始化预测器
predictor = AstroRFNet(model_path="models/net.params")

# 运行预测
labels = predictor.predict(adata, marker_genes)
result = predictor.add_predictions(adata, labels)
