#include "Tabulated.h"
#include "../../include/Exception.h"
#include "../../include/MathFunctions.h"
#include "../../include/Units.h"
#include "../../include/interp.h"
#include "../../include/interp_functions.h"
#include <cstdlib>

namespace {
std::string rns_surface_directory() {
    const char* configured = std::getenv("NS_SWORD_RNS_DIR");
    return configured && *configured ? configured : "inputs/RNS_Shapes";
}
}

TabulatedModel::TabulatedModel(double req_nounits, double mass_nounits, double omega_nounits, int model)
    : OblModelBase(req_nounits, mass_nounits, omega_nounits) {
    loadShapeFunction(model); // Load shape function data during construction
    model_name = "Tab_" + std::to_string(model);
}

TabulatedModel::TabulatedModel(double req_nounits, double mass_nounits, double omega_nounits, std::string model)
    : OblModelBase(req_nounits, mass_nounits, omega_nounits) {
    // Helper lambda to trim whitespace from a string
    loadShapeFunction(model); // Load shape function data during construction
    model_name = "Tab_" + model;
}

// Load shape function data from file and convert radii to nounits
void TabulatedModel::loadShapeFunction(int model) {
    std::string filename = rns_surface_directory() + "/" + std::to_string(model) + ".csv";
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("TabulatedModel::loadShapeFunction: Unable to open shape function file: " + filename);
    }

    double mu, radius_km;
    std::string line;
    std::getline(file, line); // Read past the header line 
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string mu_str, radius_km_str;

        // Read up to the comma delimiter
        if (std::getline(iss, mu_str, ',') && std::getline(iss, radius_km_str)) {
            try {
                double mu = std::stod(mu_str);
                double radius_km = std::stod(radius_km_str);
                // Convert radius from km to cgs and then to nounits
                double radius_nounits = Units::cgs_to_nounits(radius_km * 1.0e5, Units::LENGTH);
                
                //std::cout << mu << ", " << radius_km << ", " << radius_nounits << std::endl;
                mu_radii_map[mu] = radius_nounits;
                if (mu != 0) mu_radii_map[-mu] = radius_nounits; //Handle -1 to 0.

                //std::cout << "Parsed mu: " << mu << ", radius_nounits: " << radius_nounits << std::endl;
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid data in line: " << line << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Out of range data in line: " << line << std::endl;
            }
        } else {
            std::cerr << "Failed to parse line: " << line << std::endl;
        }
    }
    file.close();
}

// Load shape function data from file and convert radii to nounits
void TabulatedModel::loadShapeFunction(std::string model) {
    // Helper lambda to trim whitespace from a string
    auto trim = [](std::string str) {
        str.erase(0, str.find_first_not_of(" \t\n\r"));
        str.erase(str.find_last_not_of(" \t\n\r") + 1);
    };
    trim(model);

    std::string filename = rns_surface_directory() + "/" + model + ".csv";
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("TabulatedModel::loadShapeFunction: Unable to open shape function file: " + filename);
    }

    double mu, radius_km;
    std::string line;
    std::getline(file, line); // Read past the header line 
    while (std::getline(file, line)) {
        std::istringstream iss(line);
        std::string mu_str, radius_km_str;

        // Read up to the comma delimiter
        if (std::getline(iss, mu_str, ',') && std::getline(iss, radius_km_str)) {
            try {
                double mu = std::stod(mu_str);
                double radius_km = std::stod(radius_km_str);
                // Convert radius from km to cgs and then to nounits
                double radius_nounits = Units::cgs_to_nounits(radius_km * 1.0e5, Units::LENGTH);
                
                //std::cout << mu << ", " << radius_km << ", " << radius_nounits << std::endl;
                mu_radii_map[mu] = radius_nounits;
                if (mu != 0) mu_radii_map[-mu] = radius_nounits; //Handle -1 to 0.

                //std::cout << "Parsed mu: " << mu << ", radius_nounits: " << radius_nounits << std::endl;
            } catch (const std::invalid_argument& e) {
                std::cerr << "Invalid data in line: " << line << std::endl;
            } catch (const std::out_of_range& e) {
                std::cerr << "Out of range data in line: " << line << std::endl;
            }
        } else {
            std::cerr << "Failed to parse line: " << line << std::endl;
        }
    }
    file.close();
}

// TabulatedModel method using linear_interpolate
double TabulatedModel::R_at_costheta(const double& costheta) const {
    if (mu_radii_map.empty()) {
        throw std::runtime_error("TabulatedModel::R_at_costheta: No shape function data loaded.");
    }

    // Prepare data for interpolation
    std::vector<double> xp, yp;
    for (const auto& [mu, radius] : mu_radii_map) {
        xp.push_back(mu);
        yp.push_back(radius);
    }

    // Call the linear_interpolate function
    int order = 3;
    return polynomial_interpolate(xp, yp, costheta, order);
}


// Calculate dR(theta)/dtheta using: (1) SG smoothing of R(mu), (2) natural cubic derivative.
// Drop-in replacement for the previous version. Created 2025 with ChatGPT.
// Accuracy: SG (p=3) gives O(h^3) value locally; natural cubic derivative O(h^3).
// Runtime: one-time O(n * W * p^2 + n) build; per-call O(log n).
double TabulatedModel::Dtheta_R(const double& costheta) {
    // ---- Tunables (compile-time): SG window & degree ----
    // Choose odd window >= degree+1. 9/3 is a good default; increase window to smooth more.
    constexpr int SG_WINDOW = 9;
    constexpr int SG_DEGREE = 3;

    // ---- Local cache for speed across repeated calls ----
    struct Cache {
        std::vector<double> x;        // knots (mu), strictly increasing
        std::vector<double> y_raw;    // original R
        std::vector<double> y_smooth; // SG-smoothed R
        std::vector<double> M;        // second-derivatives of natural cubic on y_smooth
        std::size_t checksum = 0;
        bool built = false;
        int sg_window_used = 0;
        int sg_degree_used = 0;
    };
    static Cache cache;

    // ---- Checksum to detect any change in the table ----
    auto compute_checksum = [](const std::map<double,double>& m)->std::size_t {
        std::size_t h = 1469598103934665603ull; // FNV-1a
        auto mix = [&](std::uint64_t v){ h ^= v; h *= 1099511628211ull; };
        for (const auto& kv : m) {
            std::uint64_t a, b;
            std::memcpy(&a, &kv.first,  sizeof(double));
            std::memcpy(&b, &kv.second, sizeof(double));
            mix(a); mix(b);
        }
        mix(static_cast<std::uint64_t>(m.size()));
        return h;
    };

    const std::size_t new_sum = compute_checksum(mu_radii_map);

    // ---- Rebuild cache if table or SG params changed ----
    if (!cache.built || new_sum != cache.checksum
        || cache.sg_window_used != SG_WINDOW || cache.sg_degree_used != SG_DEGREE) {

        // Copy map → vectors (already sorted by key)
        cache.x.clear(); cache.y_raw.clear(); cache.y_smooth.clear(); cache.M.clear();
        cache.x.reserve(mu_radii_map.size());
        cache.y_raw.reserve(mu_radii_map.size());
        for (const auto& [mu, r] : mu_radii_map) {
            cache.x.push_back(mu);
            cache.y_raw.push_back(r);
        }
        const int n = static_cast<int>(cache.x.size());
        cache.y_smooth.assign(n, 0.0);
        cache.M.assign(std::max(0, n), 0.0);

        // --- Helper: adjust SG window for dataset size/degree ---
        auto adjust_window = [&](int W)->int {
            int w = W;
            if (w < 3) w = 3;
            if (w > n) w = (n % 2 ? n : n-1);
            if (w % 2 == 0) --w;
            int min_needed = SG_DEGREE + 1;
            if (w < min_needed) w = (min_needed % 2 ? min_needed : min_needed + 1);
            if (w > n) w = (n % 2 ? n : n-1);
            if (w < 3) w = 3;
            return w;
        };
        const int W = adjust_window(SG_WINDOW);
        const int p = SG_DEGREE;
        const int half = W / 2;

        // --- Savitzky–Golay smoothing on non-uniform grid (local LS) ---
        // At each knot i, fit local polynomial y ≈ sum_{k=0..p} c_k (x-x_i)^k; take smoothed value c_0.
        // (We *do not* take c_1 from SG because we want a globally C^1 result via the spline.)
        std::vector<double> c(p+1), rhs(p+1), powers(p+1);
        std::vector<double> G((p+1)*(p+1)); // normal matrix

        auto solve_small = [&](std::vector<double>& A, std::vector<double>& b, int m){ // Gaussian elim (m <= 5)
            for (int k = 0; k < m; ++k) {
                // pivot
                int piv = k;
                double best = std::fabs(A[k*m + k]);
                for (int r = k+1; r < m; ++r) {
                    double v = std::fabs(A[r*m + k]);
                    if (v > best) { best = v; piv = r; }
                }
                if (piv != k) {
                    for (int c2 = k; c2 < m; ++c2) std::swap(A[k*m + c2], A[piv*m + c2]);
                    std::swap(b[k], b[piv]);
                }
                const double diag = A[k*m + k];
                if (std::fabs(diag) < 1e-18) continue;
                // eliminate
                for (int r = k+1; r < m; ++r) {
                    const double f = A[r*m + k] / diag;
                    if (f == 0.0) continue;
                    for (int c2 = k; c2 < m; ++c2) A[r*m + c2] -= f * A[k*m + c2];
                    b[r] -= f * b[k];
                }
            }
            // back-substitute
            for (int i = m-1; i >= 0; --i) {
                double s = b[i];
                for (int j = i+1; j < m; ++j) s -= A[i*m + j]*c[j];
                const double diag = A[i*m + i];
                c[i] = (std::fabs(diag) < 1e-18) ? 0.0 : (s / diag);
            }
        };

        for (int i = 0; i < n; ++i) {
            const int i0 = std::max(0, i - half);
            const int i1 = std::min(n-1, i + half);
            const int mpts = i1 - i0 + 1;

            // Zero normal system
            std::fill(G.begin(), G.end(), 0.0);
            std::fill(rhs.begin(), rhs.end(), 0.0);

            for (int j = i0; j <= i1; ++j) {
                const double t = cache.x[j] - cache.x[i]; // center at x_i
                powers[0] = 1.0;
                for (int k = 1; k <= p; ++k) powers[k] = powers[k-1] * t;

                // Accumulate V^T V and V^T y
                for (int r = 0; r <= p; ++r) {
                    rhs[r] += powers[r] * cache.y_raw[j];
                    for (int s = 0; s <= p; ++s) {
                        G[r*(p+1) + s] += powers[r] * powers[s];
                    }
                }
            }

            // Solve G c = rhs
            std::fill(c.begin(), c.end(), 0.0);
            solve_small(G, rhs, p+1);

            cache.y_smooth[i] = c[0]; // smoothed value at x_i
        }

        // --- Build natural cubic spline on y_smooth (second derivatives M) ---
        if (n >= 3) {
            std::vector<double> a(n-2), b(n-2), cc(n-2), d(n-2);
            for (int i = 1; i <= n-2; ++i) {
                const double h_im1 = cache.x[i]   - cache.x[i-1];
                const double h_i   = cache.x[i+1] - cache.x[i];
                a[i-1] = h_im1;
                b[i-1] = 2.0 * (h_im1 + h_i);
                cc[i-1]= h_i;
                d[i-1] = 6.0 * ( (cache.y_smooth[i+1]-cache.y_smooth[i])/h_i
                               - (cache.y_smooth[i]-cache.y_smooth[i-1])/h_im1 );
            }
            // Thomas
            for (int i = 1; i <= n-3; ++i) {
                const double m = a[i] / b[i-1];
                b[i]  -= m * cc[i-1];
                d[i]  -= m * d[i-1];
            }
            std::vector<double> Min(n-2, 0.0);
            Min[n-3] = d[n-3] / b[n-3];
            for (int i = n-4; i >= 0; --i) {
                Min[i] = (d[i] - cc[i]*Min[i+1]) / b[i];
            }
            for (int i = 1; i <= n-2; ++i) cache.M[i] = Min[i-1];
            cache.M[0] = 0.0; cache.M[n-1] = 0.0;
        }

        cache.checksum = new_sum;
        cache.built = true;
        cache.sg_window_used = W;
        cache.sg_degree_used = p;
    }

    // ---- If no data, exit gracefully ----
    if (cache.x.empty()) return 0.0;

    // ---- Evaluate dR/dmu from the smoothed spline ----
    const double mu_min = cache.x.front();
    const double mu_max = cache.x.back();
    const double mu = std::min(std::max(costheta, mu_min), mu_max);

    const int n = static_cast<int>(cache.x.size());
    double dR_dmu = 0.0;

    if (n == 1) {
        dR_dmu = 0.0;
    } else if (mu <= cache.x.front()) {
        const double h = cache.x[1] - cache.x[0];
        const double A = (cache.y_smooth[1] - cache.y_smooth[0]) / h
                         - (cache.M[1] + 2.0*cache.M[0]) * h / 6.0;
        dR_dmu = A;
    } else if (mu >= cache.x.back()) {
        const double h = cache.x[n-1] - cache.x[n-2];
        const double B = (cache.y_smooth[n-1] - cache.y_smooth[n-2]) / h
                         + (2.0*cache.M[n-1] + cache.M[n-2]) * h / 6.0;
        dR_dmu = B;
    } else {
        auto it = std::upper_bound(cache.x.begin(), cache.x.end(), mu);
        const int i = static_cast<int>((it - cache.x.begin()) - 1);
        const double h  = cache.x[i+1] - cache.x[i];
        const double tL = cache.x[i+1] - mu;
        const double tR = mu - cache.x[i];
        const double Mi  = cache.M[i];
        const double Mi1 = cache.M[i+1];
        const double A = cache.y_smooth[i]   - Mi  * (h*h)/6.0;
        const double B = cache.y_smooth[i+1] - Mi1 * (h*h)/6.0;
        dR_dmu = -(Mi  * (tL*tL))/(2.0*h)
                 + (Mi1 * (tR*tR))/(2.0*h)
                 - A/h + B/h;
    }

    // ---- Chain rule to dR/dtheta ----
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu*mu));
    return -sin_theta * dR_dmu;
}


/**
// Calculate dR(theta)/dtheta using an analytic derivative from a cached natural cubic spline.
// Drop-in replacement for the previous finite-difference version.
// Created in 2025 with the help of ChatGPT.
// Convergence of derivative: O(h^3). Build O(n); per-call eval O(log n).
double TabulatedModel::Dtheta_R(const double& costheta) {
    // --- Local static cache so repeated calls are fast ---
    struct Cache {
        std::vector<double> x;   // mu knots (increasing)
        std::vector<double> y;   // R values
        std::vector<double> M;   // second derivatives at knots (natural cubic)
        std::size_t checksum = 0;
        bool built = false;
    };
    static Cache cache;

    // --- Build a simple checksum of the tabulated map to detect changes ---
    auto compute_checksum = [](const std::map<double,double>& m)->std::size_t {
        std::size_t h = 1469598103934665603ull; // FNV-1a offset basis
        auto mix = [&](std::uint64_t v){
            h ^= v; h *= 1099511628211ull;
        };
        for (const auto& kv : m) {
            std::uint64_t a, b;
            std::memcpy(&a, &kv.first,  sizeof(double));
            std::memcpy(&b, &kv.second, sizeof(double));
            mix(a); mix(b);
        }
        // also encode size in case of trailing zeros etc.
        mix(static_cast<std::uint64_t>(m.size()));
        return h;
    };

    const std::size_t new_sum = compute_checksum(mu_radii_map);

    // --- (Re)build spline if needed ---
    if (!cache.built || new_sum != cache.checksum) {
        cache.x.clear(); cache.y.clear(); cache.M.clear();

        cache.x.reserve(mu_radii_map.size());
        cache.y.reserve(mu_radii_map.size());
        for (const auto& [mu, r] : mu_radii_map) {
            cache.x.push_back(mu);
            cache.y.push_back(r);
        }

        const int n = static_cast<int>(cache.x.size());
        cache.M.assign(std::max(0, n), 0.0);

        if (n >= 3) {
            // Build natural cubic spline: solve tridiagonal for second derivatives M
            // h_{i-1} M_{i-1} + 2(h_{i-1}+h_i) M_i + h_i M_{i+1} = 6( (y_{i+1}-y_i)/h_i - (y_i - y_{i-1})/h_{i-1} )
            std::vector<double> a(n-2), b(n-2), c(n-2), d(n-2);

            for (int i = 1; i <= n-2; ++i) {
                const double h_im1 = cache.x[i]   - cache.x[i-1];
                const double h_i   = cache.x[i+1] - cache.x[i];
                a[i-1] = h_im1;
                b[i-1] = 2.0 * (h_im1 + h_i);
                c[i-1] = h_i;
                d[i-1] = 6.0 * ( (cache.y[i+1]-cache.y[i])/h_i - (cache.y[i]-cache.y[i-1])/h_im1 );
            }
            // Thomas algorithm
            for (int i = 1; i <= n-3; ++i) {
                const double m = a[i] / b[i-1];
                b[i] -= m * c[i-1];
                d[i] -= m * d[i-1];
            }
            std::vector<double> Min(n-2, 0.0);
            Min[n-3] = d[n-3] / b[n-3];
            for (int i = n-4; i >= 0; --i) {
                Min[i] = (d[i] - c[i]*Min[i+1]) / b[i];
            }
            for (int i = 1; i <= n-2; ++i) cache.M[i] = Min[i-1];
            cache.M[0] = 0.0; cache.M[n-1] = 0.0; // natural boundaries
        }
        // n==0/1/2 → M already zeroed; linear/constant cases handled below.

        cache.checksum = new_sum;
        cache.built = true;
    }

    // --- Handle empty data gracefully ---
    if (cache.x.empty()) return 0.0;

    // Clamp mu to data domain to avoid extrapolation noise
    const double mu_min = cache.x.front();
    const double mu_max = cache.x.back();
    const double mu = std::min(std::max(costheta, mu_min), mu_max);

    // Binary search for interval [x[i], x[i+1]]
    const int n = static_cast<int>(cache.x.size());
    double dR_dmu = 0.0;

    if (n == 1) {
        dR_dmu = 0.0;
    } else if (mu <= cache.x.front()) {
        // Left endpoint derivative from first interval
        const double h = cache.x[1] - cache.x[0];
        const double A = (cache.y[1] - cache.y[0]) / h - (cache.M[1] + 2.0*cache.M[0]) * h / 6.0;
        dR_dmu = A;
    } else if (mu >= cache.x.back()) {
        // Right endpoint derivative from last interval
        const double h = cache.x[n-1] - cache.x[n-2];
        const double B = (cache.y[n-1] - cache.y[n-2]) / h + (2.0*cache.M[n-1] + cache.M[n-2]) * h / 6.0;
        dR_dmu = B;
    } else {
        auto it = std::upper_bound(cache.x.begin(), cache.x.end(), mu);
        const int i = static_cast<int>((it - cache.x.begin()) - 1);
        const double h  = cache.x[i+1] - cache.x[i];
        const double tL = cache.x[i+1] - mu;   // left distance
        const double tR = mu - cache.x[i];     // right distance
        const double Mi  = cache.M[i];
        const double Mi1 = cache.M[i+1];
        const double A = cache.y[i]   - Mi  * (h*h)/6.0;
        const double B = cache.y[i+1] - Mi1 * (h*h)/6.0;

        // Analytic derivative of natural cubic on [x_i, x_{i+1}]:
        // S'(mu) = -Mi*(tL^2)/(2h) + Mi1*(tR^2)/(2h) - A/h + B/h
        dR_dmu = -(Mi  * (tL*tL))/(2.0*h)
                 + (Mi1 * (tR*tR))/(2.0*h)
                 - A/h + B/h;
    }

    // Chain rule: mu = cos(theta), dmu/dtheta = -sin(theta) = -sqrt(1 - mu^2) on θ∈[0,π]
    const double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu*mu));
    return -sin_theta * dR_dmu;
}
*/


// External function to create the model
extern "C" OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model) {
    return new TabulatedModel(req_nounits, mass_nounits, omega_nounits, model);
}

extern "C" OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name) {
    return new TabulatedModel(req_nounits, mass_nounits, omega_nounits, model_name);
}
