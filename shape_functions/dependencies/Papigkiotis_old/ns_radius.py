import sys
import numpy as np
import torch
from pathlib import Path
import sys

CANONICAL_DIR = Path(__file__).resolve().parents[1] / "Papigkiotis"
sys.path.insert(0, str(CANONICAL_DIR))
from DNN import Regressor

abs_mu_min, abs_mu_max = 0., 1. 
C_min, C_max = 0.0876346858172578, 0.3094541325480277
sigma_min, sigma_max = 0., 0.9612274013913829
eccentricity_min, eccentricity_max = 0., 0.7797886226038347

min_values = np.array([abs_mu_min, C_min, sigma_min, eccentricity_min])
max_values = np.array([abs_mu_max, C_max, sigma_max, eccentricity_max])
feature_scaler = lambda data: (data - min_values) / (max_values - min_values)
#no_scaler = lambda data: data

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# ---------- load ANN --------------------------------------------------------
regressor_shape = Regressor(input_dimension=4, feature_scaler=feature_scaler).to(device) #Regressor(input_dimension=4, feature_scaler=no_scaler).to(device) #
weights_path = Path("shape_functions/dependencies/Papigkiotis/Model/Surface/Surface-model.pth")
regressor_shape.load_state_dict(torch.load(weights_path, map_location=device))
regressor_shape.eval()

# ---------- analytic helpers ------------------------------------------------
def r_mu(C, sigma, R_pole, R_eq, num = 20):
    mu = np.linspace(0, 1, num=num, dtype=np.float32)
    C_np = np.array([C for _ in range(0, num)], dtype=np.float32)
    sigma_np = np.array([sigma for _ in range(0, num)], dtype=np.float32)
    e_np = np.array([np.sqrt(1 - np.square(R_pole / R_eq)) for _ in range(0, num)], dtype=np.float32)
    data = np.array([mu, C_np, sigma_np, e_np])
    x = torch.tensor(data).T.to(device)
    model_estimation = regressor_shape(x) * (R_eq - R_pole) + R_pole
    model_estimation = model_estimation.cpu().detach().numpy()
    
    model_estimation = model_estimation.ravel().astype(np.float32)
    
    return mu, model_estimation

def R_pole(R_e, C, sigma):
    """
    Compute the polar radius R_pole given equatorial radius R_e, compactness C, and reduced spin sigma,
    using the coefficients from Table IV (arXiv:2501.18544 Eq. 25).
    """
    # Coefficients from Table IV
    A = {
        (0, 0):  0.942328,
        (0, 1): -0.617711,
        (0, 2):  0.544639,
        (0, 3): -0.440968,
        (0, 4):  0.196118,
        (1, 0):  1.296632,
        (1, 1): -1.458921,
        (1, 2): -0.226904,
        (1, 3):  0.527775,
        (2, 0): -10.45611,
        (2, 1):  8.668382,
        (2, 2): -2.506686,
        (3, 0):  36.131881,
        (3, 1): -7.524662,
        (4, 0): -45.301523,
    }
    # Evaluate the polynomial
    R_pole_over_R_e = 0.0
    for n in range(5):
        for m in range(5-n):
            coeff = A.get((n, m), 0.0)
            R_pole_over_R_e += coeff * (C**n) * (sigma**m)
    return R_e * R_pole_over_R_e

def R_pole_notebook_direct(R_e, C, sigma):
    return (
        R_e*(68.7151320054541 * C**5 -117.228807812831 * C**4 * sigma -
        71.6699195745365 * C**4 - 8.83436790963497 * C**3 * sigma**2 +
        90.4556842423828 * C**3 * sigma + 28.980101203154 * C**3 -
        2.1673404784385 * C**2 * sigma**3 + 4.94947873047332 * C**2 * sigma**2 -
        21.6125323056704 * C**2 * sigma - 5.6347743027431 * C**2 + 0.274033536311414 * C * sigma**4 +
        1.01941492279495 * C * sigma**3 - 1.95940441780677 * C * sigma**2 +
        2.60894540594265 * C * sigma + 0.522844117974324 * C + 0.627942468515494 * sigma**5 -
        1.16678579309014 * sigma**4 + 0.508122852358802 * sigma**3 + 0.371402374708614 * sigma**2 -
        0.784944132327725 * sigma + 0.981279105367004)
    )

# ---------- CLI -------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 4:
        raise SystemExit(
            "usage: python ns_radius.py R_eq C sigma"
        )
    
    R_eq = float(sys.argv[1]) #0.7
    Compactness = float(sys.argv[2])
    sigma = float(sys.argv[3]) #0.8
    
    R_p = R_pole(R_eq, Compactness, sigma)
    #R_p_notebook = R_pole_notebook_direct(R_eq, Compactness, sigma)
    #print(Compactness, sigma, R_p, R_p_notebook, R_eq)
    mu, R_res = r_mu(Compactness, sigma, R_p, R_eq)
    
    #print(mu[-1], ",", R_res[-1])

    #print((R_res[-1]-R_p)/R_p*100)

    for (mu_i, R_i) in zip(mu, R_res):
        print(mu_i, ",", R_i)
    
    #print(R_pole)
