import torch
import torch.nn as nn
import torch.optim as optim

class Regressor(nn.Module):
    def __init__(self, input_dimension, feature_scaler):
        super().__init__()
        negative_slope = 0.1
        self.feature_scaler = feature_scaler

        self.MLP = nn.Sequential(
            nn.Linear(input_dimension, 200),
            # nn.BatchNorm1d(200),
            nn.LeakyReLU(negative_slope=negative_slope, inplace=True),

            nn.Linear(200, 100),
            # nn.BatchNorm1d(100),
            nn.LeakyReLU(negative_slope=negative_slope, inplace=True),

            nn.Linear(100, 50),
            # nn.BatchNorm1d(50),
            nn.LeakyReLU(negative_slope=negative_slope, inplace=True),

            nn.Linear(50, 25),
            # nn.BatchNorm1d(25),
            nn.LeakyReLU(negative_slope=negative_slope, inplace=True),

            nn.Linear(25, 10),
            # nn.BatchNorm1d(25),
            nn.LeakyReLU(negative_slope=negative_slope, inplace=True),

            nn.Linear(10, 1),
            nn.Sigmoid(),
            )

    def forward(self, x):
        x = self.MLP(x)
        return x

    def predict(self, x):
        x = self.feature_scaler(x).astype('float32')
        x = torch.from_numpy(x).to(self.device)
        x = self.MLP(x)
        x = x.cpu().detach().numpy()
        return x

    def predict_tr(self, x):
        return self.MLP(x)

    def set_device(self, device):
        self.device = device
