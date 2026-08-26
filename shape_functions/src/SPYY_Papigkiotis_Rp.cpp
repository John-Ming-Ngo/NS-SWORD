#include "SPYY_Papigkiotis_Rp.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <vector>

namespace {
const char* kCoefficients = "inputs/shape_functions/Silva_Eliptical_Coefficients.csv";

double papigkiotis_polar_ratio(double C, double sigma) {
  const double A[5][5] = {
      {0.942328, -0.617711, 0.544639, -0.440968, 0.196118},
      {1.296632, -1.458921, -0.226904, 0.527775, 0.0},
      {-10.45611, 8.668382, -2.506686, 0.0, 0.0},
      {36.131881, -7.524662, 0.0, 0.0, 0.0},
      {-45.301523, 0.0, 0.0, 0.0, 0.0}};
  double ratio = 0.0;
  for (int n = 0; n <= 4; ++n)
    for (int m = 0; m <= 4 - n; ++m)
      ratio += A[n][m] * std::pow(C, n) * std::pow(sigma, m);
  return ratio;
}
}

double SPYYPapigkiotisRp::compactness(double mass, double req) { return mass / req; }
double SPYYPapigkiotisRp::spin(double omega, double mass, double req) {
  return omega * omega * req * req * req / mass;
}

SPYYPapigkiotisRp::SPYYPapigkiotisRp(double req, double mass, double omega, int)
    : OblModelBase(req, mass, omega), C_(compactness(mass, req)),
      sigma_(spin(omega, mass, req)), polar_ratio_(papigkiotis_polar_ratio(C_, sigma_)),
      e2_(0.0), a2_(0.0), a4_(0.0) {
  model_name = "SPYY_Papigkiotis_Rp";
  if (!(req > 0.0) || !(mass > 0.0) || !std::isfinite(C_) ||
      !(sigma_ >= 0.0) || !std::isfinite(sigma_))
    throw std::domain_error("SPYY_Papigkiotis_Rp: invalid Re, M, or spin.");
  if (!(polar_ratio_ > 0.0) || !std::isfinite(polar_ratio_))
    throw std::domain_error("SPYY_Papigkiotis_Rp: fitted Rp/Re is nonpositive or non-finite.");
  e2_ = std::max(0.0, 1.0 - polar_ratio_ * polar_ratio_);
  load_spyy_coefficients();
  a2_ = evaluate_coefficient(a2_coeff_);
  a4_ = evaluate_coefficient(a4_coeff_);
  if (!std::isfinite(a2_) || !std::isfinite(a4_))
    throw std::domain_error("SPYY_Papigkiotis_Rp: non-finite SPYY coefficient.");
}

SPYYPapigkiotisRp::SPYYPapigkiotisRp(double req, double mass, double omega, std::string)
    : SPYYPapigkiotisRp(req, mass, omega, 0) {}

void SPYYPapigkiotisRp::load_spyy_coefficients() {
  const int wanted = sigma_ <= 0.2 ? 1 : 2;
  std::ifstream file(kCoefficients);
  if (!file) throw std::runtime_error(std::string("Unable to open file: ") + kCoefficients);
  std::string line;
  std::getline(file, line);
  while (std::getline(file, line)) {
    std::stringstream stream(line);
    std::vector<std::string> token;
    std::string item;
    while (std::getline(stream, item, ',')) token.push_back(item);
    if (token.size() < 23) throw std::runtime_error("Malformed SPYY coefficient row.");
    if (std::stoi(token[0]) != wanted) continue;
    for (int i = 0; i < 7; ++i) {
      a2_coeff_[i] = std::stod(token[9 + i]);
      a4_coeff_[i] = std::stod(token[16 + i]);
    }
    return;
  }
  throw std::runtime_error("Required SPYY coefficient row not found.");
}

double SPYYPapigkiotisRp::evaluate_coefficient(const double c[7]) const {
  return c[0] + c[1] * std::sqrt(sigma_) + c[2] * sigma_ + c[3] * C_ +
         c[4] * sigma_ * C_ + c[5] * sigma_ * sigma_ + c[6] * C_ * C_;
}

double SPYYPapigkiotisRp::g(double mu) const {
  const double mu2 = mu * mu;
  return 1.0 + a2_ * mu2 + a4_ * mu2 * mu2 - (1.0 + a2_ + a4_) * mu2 * mu2 * mu2;
}

double SPYYPapigkiotisRp::dg_dmu(double mu) const {
  return 2.0 * a2_ * mu + 4.0 * a4_ * mu * mu * mu -
         6.0 * (1.0 + a2_ + a4_) * std::pow(mu, 5);
}

double SPYYPapigkiotisRp::validated_radical(double mu) const {
  if (!std::isfinite(mu) || std::fabs(mu) > 1.0 + 1e-12)
    throw std::domain_error("SPYY_Papigkiotis_Rp: cos(theta) outside [-1, 1].");
  mu = std::max(-1.0, std::min(1.0, mu));
  const double denominator = 1.0 - e2_ * g(mu);
  if (!(denominator > 0.0) || !std::isfinite(denominator))
    throw std::domain_error("SPYY_Papigkiotis_Rp: invalid surface denominator.");
  const double argument = (1.0 - e2_) / denominator;
  if (!(argument >= 0.0) || !std::isfinite(argument))
    throw std::domain_error("SPYY_Papigkiotis_Rp: invalid square-root argument.");
  return argument;
}

double SPYYPapigkiotisRp::R_at_costheta(const double& mu) const {
  return Req_nounits() * std::sqrt(validated_radical(mu));
}

double SPYYPapigkiotisRp::Dtheta_R(const double& input_mu) {
  const double argument = validated_radical(input_mu);
  const double mu = std::max(-1.0, std::min(1.0, input_mu));
  const double denominator = 1.0 - e2_ * g(mu);
  const double sin_theta = std::sqrt(std::max(0.0, 1.0 - mu * mu));
  const double radius = Req_nounits() * std::sqrt(argument);
  return -sin_theta * radius * e2_ * dg_dmu(mu) / (2.0 * denominator);
}

double SPYYPapigkiotisRp::polar_ratio() const { return polar_ratio_; }

extern "C" OblModelBase* createOblModelBase(const double& req, const double& mass,
                                              const double& omega, int model) {
  return new SPYYPapigkiotisRp(req, mass, omega, model);
}

extern "C" OblModelBase* createOblModelBaseStr(const double& req, const double& mass,
                                                 const double& omega, std::string model) {
  return new SPYYPapigkiotisRp(req, mass, omega, model);
}
