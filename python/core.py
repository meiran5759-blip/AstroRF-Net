import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
from scipy.sparse import issparse
from skorch import NeuralNetClassifier
import anndata

class AstroRFNet:
    """AstroRF-Net 细胞类型预测算法的完整实现"""
    
    @staticmethod
    def convert_to_cpm(adata):
        """Convert raw counts to log2(CPM+1)"""
        counts_matrix = adata.X
        if issparse(counts_matrix):
            counts_matrix = counts_matrix.toarray()
        total_counts = np.maximum(counts_matrix.sum(axis=1), 1)
        cpm = (counts_matrix / total_counts[:, None]) * 1e6
        return anndata.AnnData(np.log2(cpm + 1), obs=adata.obs, var=adata.var, dtype='float32')

    @staticmethod
    def prepare_tensor(adata, marker_genes):
        """Prepare PyTorch tensor from marker genes"""
        valid_genes = sorted(g for g in marker_genes if g in adata.var_names)
        adata_filtered = adata[:, valid_genes].copy()
        X = adata_filtered.X.toarray() if issparse(adata_filtered.X) else adata_filtered.X
        return torch.from_numpy(X).float(), valid_genes

    class NeuralNet(nn.Module):
        """神经网络架构定义"""
        def __init__(self, input_size=50, hidden_sizes=(400, 200, 50), output_size=3):
            super().__init__()
            self.layers = nn.ModuleList()
            prev_size = input_size
            for size in hidden_sizes:
                self.layers.append(nn.Linear(prev_size, size))
                prev_size = size
            self.output = nn.Linear(prev_size, output_size)
            
        def forward(self, x):
            for layer in self.layers:
                x = F.relu(layer(x))
            return F.softmax(self.output(x), dim=1)

    def __init__(self, model_path, device=None):
        """初始化预训练模型"""
        self.device = device or ('cuda' if torch.cuda.is_available() else 'cpu')
        self.model = self._load_model(model_path)
        
    def _load_model(self, path):
        """加载预训练模型"""
        model = NeuralNetClassifier(
            self.NeuralNet,
            device=self.device,
            optimizer=torch.optim.SGD,
            lr=0.1,
            max_epochs=1800,
            batch_size=500
        )
        model.initialize()
        model.module_.load_state_dict(torch.load(path, map_location=self.device))
        return model

    def predict(self, adata, marker_genes):
        """完整预测流程"""
        # 数据预处理
        adata_cpm = self.convert_to_cpm(adata)
        X_tensor, _ = self.prepare_tensor(adata_cpm, marker_genes)
        
        # 预测
        with torch.no_grad():
            y_pred = self.model.predict(X_tensor)
        
        return y_pred

    @staticmethod
    def add_predictions(adata, labels, col_name="AstroRF-Net"):
        """将预测结果添加到anndata对象"""
        celltype_map = {0: "Astrocyte", 1: "Glioblast", 2: "Radial-Glial"}
        adata.obs[col_name] = [celltype_map[x] for x in labels]
        return adata
