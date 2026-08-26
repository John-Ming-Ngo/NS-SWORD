"""PVLS surface-radius entry point with feature scaling enabled."""

import sys
import importlib.util
from pathlib import Path

import numpy as np
import torch

LEGACY_DIR = Path(__file__).resolve().parents[1] / "Papigkiotis_old"
sys.path.insert(0, str(LEGACY_DIR))
_LEGACY_SPEC = importlib.util.spec_from_file_location(
    "papigkiotis_legacy_ns_radius", LEGACY_DIR / "ns_radius.py"
)
legacy = importlib.util.module_from_spec(_LEGACY_SPEC)
_LEGACY_SPEC.loader.exec_module(legacy)


def r_mu(C, sigma, R_pole, R_eq, num=20):
    mu = np.linspace(0.0, 1.0, num=num, dtype=np.float32)
    eccentricity = np.sqrt(np.clip(1.0 - np.square(R_pole / R_eq), 0.0, None))
    data = np.column_stack(
        [
            mu,
            np.full_like(mu, C),
            np.full_like(mu, sigma),
            np.full_like(mu, eccentricity),
        ]
    )
    scaled = legacy.feature_scaler(data).astype(np.float32)
    x = torch.tensor(scaled, dtype=torch.float32, device=legacy.device)
    with torch.no_grad():
        estimate = legacy.regressor_shape(x).cpu().numpy().ravel()
    # The network approximates a normalized profile, but its endpoint values
    # are not mathematically exact.  Affinely normalize the whole prediction
    # so the supplied equatorial and polar radii remain the actual endpoints;
    # this removes the smooth, model-wide radius offset without inserting a
    # synthetic grid point.
    denominator = float(estimate[0] - estimate[-1])
    if not np.isfinite(denominator) or abs(denominator) < np.finfo(float).eps:
        raise ValueError("PVLS surface prediction has degenerate endpoints")
    normalized = (estimate.astype(np.float64) - float(estimate[-1])) / denominator
    radius = normalized * (R_eq - R_pole) + R_pole
    radius[0] = R_eq
    radius[-1] = R_pole
    return mu, radius


if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise SystemExit("usage: python ns_radius.py R_eq C sigma")
    R_eq, compactness, sigma = map(float, sys.argv[1:4])
    R_pole = legacy.R_pole(R_eq, compactness, sigma)
    mu, radii = r_mu(compactness, sigma, R_pole, R_eq)
    for mu_i, radius_i in zip(mu, radii):
        print(mu_i, ",", radius_i)
