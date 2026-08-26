#ifndef PVLS_CACHE_H
#define PVLS_CACHE_H

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#ifndef PAPIGKIOTIS_ENFORCE_TRAINED_DOMAIN
#define PAPIGKIOTIS_ENFORCE_TRAINED_DOMAIN 1
#endif

#ifndef PAPIGKIOTIS_ENABLE_CACHE
#define PAPIGKIOTIS_ENABLE_CACHE PAPIGKIOTIS_ENFORCE_TRAINED_DOMAIN
#endif

namespace PVLSCache {

static const double C_MIN = 0.0876346858172578;
static const double C_MAX = 0.3094541325480277;
static const double SIGMA_MIN = 0.0;
static const double SIGMA_MAX = 0.9612274013913829;
static const double ECCENTRICITY_MIN = 0.0;
static const double ECCENTRICITY_MAX = 0.7797886226038347;
static const double RP_REQ_MIN = 0.6260428931452016;
static const double RP_REQ_MAX = 1.0;

inline double polar_ratio(double compactness, double sigma) {
  static const double coefficients[5][5] = {
      { 0.942328, -0.617711,  0.544639, -0.440968,  0.196118},
      { 1.296632, -1.458921, -0.226904,  0.527775,  0.0},
      {-10.45611,  8.668382, -2.506686,  0.0,       0.0},
      {36.131881, -7.524662,  0.0,       0.0,       0.0},
      {-45.301523, 0.0,       0.0,       0.0,       0.0}
  };
  double ratio = 0.0;
  for (int n = 0; n < 5; ++n) {
    for (int m = 0; m < 5 - n; ++m) {
      ratio += coefficients[n][m] * std::pow(compactness, n) * std::pow(sigma, m);
    }
  }
  return ratio;
}

inline bool in_trained_domain(double compactness, double sigma,
                              double* ratio_out = NULL,
                              double* eccentricity_out = NULL) {
  const double ratio = polar_ratio(compactness, sigma);
  const double eccentricity_squared = 1.0 - ratio * ratio;
  const double eccentricity = eccentricity_squared >= 0.0
      ? std::sqrt(eccentricity_squared)
      : std::numeric_limits<double>::quiet_NaN();
  if (ratio_out) *ratio_out = ratio;
  if (eccentricity_out) *eccentricity_out = eccentricity;
  return std::isfinite(compactness) && std::isfinite(sigma) &&
         std::isfinite(ratio) && std::isfinite(eccentricity) &&
         compactness >= C_MIN && compactness <= C_MAX &&
         sigma >= SIGMA_MIN && sigma <= SIGMA_MAX &&
         eccentricity >= ECCENTRICITY_MIN && eccentricity <= ECCENTRICITY_MAX &&
         ratio >= RP_REQ_MIN && ratio <= RP_REQ_MAX;
}

inline void require_trained_domain(double compactness, double sigma) {
#if PAPIGKIOTIS_ENFORCE_TRAINED_DOMAIN
  double ratio = 0.0;
  double eccentricity = 0.0;
  if (!in_trained_domain(compactness, sigma, &ratio, &eccentricity)) {
    std::ostringstream message;
    message << "PVLS: requested model is outside the complete trained feature domain"
            << " (C=" << compactness << ", sigma=" << sigma
            << ", Rp/Req=" << ratio << ", eccentricity=" << eccentricity << ").";
    throw std::domain_error(message.str());
  }
#else
  (void)compactness;
  (void)sigma;
#endif
}

#if PAPIGKIOTIS_ENABLE_CACHE

static const char* const MAGIC = "PVLS_INFERENCE_CACHE_V1";
static const std::size_t SURFACE_COUNT = 20;
static const std::size_t DERIVATIVE_COUNT = 500;

inline std::string rendered_key(double req, double compactness, double sigma) {
  // Match the full-precision values passed to the Python fallback command.
  std::ostringstream key;
  key << std::setprecision(17) << req << '\t' << compactness << '\t' << sigma;
  return key.str();
}

inline std::uint64_t fnv1a_update(std::uint64_t hash, const char* data, std::size_t size) {
  const std::uint64_t prime = UINT64_C(1099511628211);
  for (std::size_t i = 0; i < size; ++i) {
    hash ^= static_cast<unsigned char>(data[i]);
    hash *= prime;
  }
  return hash;
}

inline std::string hex_hash(std::uint64_t hash) {
  std::ostringstream rendered;
  rendered << std::hex << std::setfill('0') << std::setw(16) << hash;
  return rendered.str();
}

inline std::string key_hash(const std::string& value) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  return hex_hash(fnv1a_update(hash, value.data(), value.size()));
}

inline std::string dependency_fingerprint() {
  static const char* const dependencies[] = {
      "shape_functions/dependencies/Papigkiotis/ns_radius.py",
      "shape_functions/dependencies/Papigkiotis/ns_log_derivative.py",
      "shape_functions/dependencies/Papigkiotis_old/ns_radius.py",
      "shape_functions/dependencies/Papigkiotis_old/ns_log_derivative.py",
      "shape_functions/dependencies/Papigkiotis/DNN.py",
      "shape_functions/dependencies/Papigkiotis/Model/Surface/Surface-model.pth",
      "shape_functions/dependencies/Papigkiotis/Model/Derivative/Derivative-model.pth"
  };
  std::uint64_t hash = UINT64_C(14695981039346656037);
  char buffer[8192];
  for (std::size_t i = 0; i < sizeof(dependencies) / sizeof(dependencies[0]); ++i) {
    const std::string path(dependencies[i]);
    hash = fnv1a_update(hash, path.data(), path.size());
    const char separator = '\0';
    hash = fnv1a_update(hash, &separator, 1);
    std::ifstream input(path.c_str(), std::ios::in | std::ios::binary);
    if (!input) throw std::runtime_error("PVLS cache: cannot fingerprint dependency " + path);
    while (input) {
      input.read(buffer, sizeof(buffer));
      const std::streamsize count = input.gcount();
      if (count > 0) hash = fnv1a_update(hash, buffer, static_cast<std::size_t>(count));
    }
  }
  return hex_hash(hash);
}

inline std::string path_for_key(const std::string& key) {
  const char* configured = std::getenv("PVLS_CACHE_DIR");
  const std::string directory = configured && *configured
      ? std::string(configured) : std::string("build_objs/cache/PVLS");
  const char last = directory.empty() ? '\0' : directory[directory.size() - 1];
  return directory + ((last == '/' || last == '\\') ? "" : "/") + key_hash(key) + ".pvls";
}

inline void mirror(std::map<double, double>* values, bool antisymmetric) {
  if (values->empty() || values->begin()->first <= -1.0) return;
  for (std::map<double, double>::const_iterator it = values->begin(); it != values->end(); ++it) {
    if (it->first > 0.0) {
      values->emplace(-it->first, antisymmetric ? -it->second : it->second);
    }
  }
}

inline bool load(double req, double compactness, double sigma,
                 std::map<double, double>* surface,
                 std::map<double, double>* derivative) {
  const std::string key = rendered_key(req, compactness, sigma);
  const std::string path = path_for_key(key);
  std::ifstream input(path.c_str());
  if (!input) return false;

  std::string magic, stored_key, domain_marker, fingerprint_marker;
  if (!std::getline(input, magic) || !std::getline(input, stored_key) ||
      !std::getline(input, domain_marker) || !std::getline(input, fingerprint_marker)) {
    throw std::runtime_error("PVLS cache: truncated header in " + path);
  }
  if (magic != MAGIC || stored_key != key || domain_marker != "trained-domain") {
    throw std::runtime_error("PVLS cache: identity or domain validation failed in " + path);
  }
  if (fingerprint_marker != "fingerprint\t" + dependency_fingerprint()) {
    throw std::runtime_error("PVLS cache: inference dependency fingerprint mismatch in " + path);
  }

  std::size_t surface_count = 0, derivative_count = 0;
  if (!(input >> surface_count >> derivative_count) ||
      surface_count != SURFACE_COUNT || derivative_count != DERIVATIVE_COUNT) {
    throw std::runtime_error("PVLS cache: unexpected table dimensions in " + path);
  }
  surface->clear();
  derivative->clear();
  for (std::size_t i = 0; i < surface_count + derivative_count; ++i) {
    char table = '\0';
    double mu = 0.0, value = 0.0;
    if (!(input >> table >> mu >> value) || !std::isfinite(mu) || !std::isfinite(value) ||
        mu < 0.0 || mu > 1.0) {
      throw std::runtime_error("PVLS cache: malformed inference row in " + path);
    }
    std::map<double, double>* destination = table == 'S' ? surface
        : (table == 'D' ? derivative : static_cast<std::map<double, double>*>(NULL));
    if (!destination || !destination->emplace(mu, value).second) {
      throw std::runtime_error("PVLS cache: duplicate or invalid inference row in " + path);
    }
  }
  if (surface->size() != surface_count || derivative->size() != derivative_count) {
    throw std::runtime_error("PVLS cache: incomplete inference tables in " + path);
  }
  std::string trailing;
  if (input >> trailing) {
    throw std::runtime_error("PVLS cache: unexpected trailing content in " + path);
  }
  mirror(surface, false);
  mirror(derivative, true);
  return true;
}

#else

inline bool load(double, double, double,
                 std::map<double, double>*, std::map<double, double>*) {
  return false;
}

#endif // PAPIGKIOTIS_ENABLE_CACHE

} // namespace PVLSCache

#endif // PVLS_CACHE_H
