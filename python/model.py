import os
import torch
import torch.nn as nn
import torch.nn.functional as F
from skorch import NeuralNetClassifier

DEVICE_DEFAULT = "cuda" if torch.cuda.is_available() else "cpu"

class NeuralNet(nn.Module):
    def __init__(self, input_size, hidden_size1, hidden_size2, hidden_size3, output_size):
        super().__init__()
        self.fc1 = nn.Linear(input_size, hidden_size1)
        self.fc2 = nn.Linear(hidden_size1, hidden_size2)
        self.fc3 = nn.Linear(hidden_size2, hidden_size3)
        self.fc4 = nn.Linear(hidden_size3, output_size)

    def forward(self, x):
        x = F.relu(self.fc1(x))
        x = F.relu(self.fc2(x))
        x = F.relu(self.fc3(x))
        x = self.fc4(x)
        return F.softmax(x, dim=1)

def build_model(
    device: str = DEVICE_DEFAULT,
    input_size: int = 50,
    h1: int = 400,
    h2: int = 200,
    h3: int = 50,
    output_size: int = 3,
    lr: float = 0.1
) -> NeuralNetClassifier:
    model = NeuralNetClassifier(
        NeuralNet,
        iterator_train__shuffle=True,
        device=device,
        optimizer=torch.optim.SGD,
        lr=lr,
        max_epochs=1800,
        module__input_size=input_size,
        module__hidden_size1=h1,
        module__hidden_size2=h2,
        module__hidden_size3=h3,
        module__output_size=output_size,
        batch_size=500
    )
    model.initialize()
    return model

def load_model_weights(model: NeuralNetClassifier, param_path: str, device: str = DEVICE_DEFAULT):
    if not os.path.exists(param_path):
        raise FileNotFoundError(f"Model weights not found: {param_path}")
    state = torch.load(param_path, map_location=device)
    model.module_.load_state_dict(state)
    model.module_.eval()
    model.module_.to(device)
    return model
