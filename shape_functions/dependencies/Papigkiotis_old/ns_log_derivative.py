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
Rpol_Req_min, Rpol_Req_max = 0.6260428931452016, 1.0

min_values = np.array([abs_mu_min, C_min, sigma_min, Rpol_Req_min])
max_values = np.array([abs_mu_max, C_max, sigma_max, Rpol_Req_max])

feature_scaler = lambda data: (data - min_values) / (max_values - min_values)

device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

# ---------- load ANN --------------------------------------------------------
model_path = 'shape_functions/dependencies/Papigkiotis/Model/Derivative/Derivative-model.pth'

regressor = Regressor(input_dimension=4, feature_scaler=feature_scaler).to(device)
regressor.set_device(device)
regressor.load_state_dict(torch.load(model_path, map_location=torch.device(device)))
regressor.eval()


# ---------- analytic helpers ---------------------------------
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

def dLogR_dtheta_Max(C, sigma, R):
    return (
        -0.403276601295374 * C**3 +0.661519097793562 * C**2 * sigma +
        2.57434191155828 * C**2 * R -2.23108454811518 * C**2 -
        8.67625739742325 * C * sigma**2 - 39.6796750693722 * C * sigma * R +
        38.4663786235779 * C * sigma - 45.0567823427512 * C * R**2 +
        87.0956609308423 * C * R - 42.1349805655397 * C +
        4.07991541887412 * sigma**3 + 30.4880078788823 * sigma**2 * R -
        27.2418077134724 * sigma**2 + 72.1353182129003 * sigma * R**2 -
        130.455828718591 * sigma * R + 58.5590616356724 * sigma +
        53.9663039265715 * R**3 - 146.858522848505 * R**2 +
        131.341936932855 * R - 38.4414384733067
    )

def dlogR_dmu(C, sigma, R_pole, R_eq, num=500):
    """
    Returns (mu, dlnR/dθ) on a cosine-latitude grid.
    """
    mu   = np.linspace(0.0, 1.0, num, dtype=np.float32)
    C_np = np.full_like(mu, C)
    sig_np = np.full_like(mu, sigma)
    eccentricity = np.array([np.sqrt(1 - np.square(R_pole / R_eq)) for _ in range(0, num)], dtype=np.float32)
    RpReq = np.full_like(mu, R_pole / R_eq)

    x = torch.tensor(np.column_stack([mu, C_np, sig_np, RpReq]),
                     dtype=torch.float32, device=device)

    with torch.no_grad():
        y_hat = regressor(x).cpu().numpy().ravel().astype(np.float32)

    dmax = dLogR_dtheta_Max(C, sigma, R_p/R_eq)
    #print(dmax)
    dlogRdmu = y_hat * dmax     # ANN output ∈ (0,1)
    return mu, dlogRdmu

# ---------- CLI -------------------------------------------------------------
if __name__ == "__main__":
    if len(sys.argv) < 4:
        sys.exit("usage: python ns_log_derivative.py  R_eq  C  sigma")

    R_eq  = float(sys.argv[1]) #0.7 
    C     = float(sys.argv[2]) #0.15
    sigma = float(sys.argv[3]) #0.8

    R_p = R_pole(R_eq, C, sigma)
    mu, dlnR = dlogR_dmu(C, sigma, R_p, R_eq)

    for m, d in zip(mu, dlnR):
        print(f"{m} , {d}")
        
