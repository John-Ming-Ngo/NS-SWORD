"""PVLS log-radius derivative entry point with feature scaling enabled."""

import sys
import importlib.util
from pathlib import Path

import numpy as np
import torch

LEGACY_DIR = Path(__file__).resolve().parents[1] / "Papigkiotis_old"
sys.path.insert(0, str(LEGACY_DIR))
_LEGACY_SPEC = importlib.util.spec_from_file_location(
    "papigkiotis_legacy_ns_log_derivative", LEGACY_DIR / "ns_log_derivative.py"
)
legacy = importlib.util.module_from_spec(_LEGACY_SPEC)
_LEGACY_SPEC.loader.exec_module(legacy)


def dlogR_dmu(C, sigma, R_pole, R_eq, num=500):
    mu = np.linspace(0.0, 1.0, num=num, dtype=np.float32)
    data = np.column_stack(
        [
            mu,
            np.full_like(mu, C),
            np.full_like(mu, sigma),
            np.full_like(mu, R_pole / R_eq),
        ]
    )
    scaled = legacy.feature_scaler(data).astype(np.float32)
    x = torch.tensor(scaled, dtype=torch.float32, device=legacy.device)
    with torch.no_grad():
        estimate = legacy.regressor(x).cpu().numpy().ravel().astype(np.float32)
    maximum = legacy.dLogR_dtheta_Max(C, sigma, R_pole / R_eq)
    return mu, estimate * maximum


if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise SystemExit("usage: python ns_log_derivative.py R_eq C sigma")
    R_eq, compactness, sigma = map(float, sys.argv[1:4])
    R_pole = legacy.R_pole(R_eq, compactness, sigma)
    mu, derivatives = dlogR_dmu(compactness, sigma, R_pole, R_eq)
    for mu_i, derivative_i in zip(mu, derivatives):
        print(f"{mu_i} , {derivative_i}")
