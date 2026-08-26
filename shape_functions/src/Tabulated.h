#ifndef TABULATED_H
#define TABULATED_H

#include "../../include/OblModelBase.h"
#include <exception>
#include <string>
#include <map>
#include <cstdint>    // for std::uint64_t
#include <algorithm>  // for std::upper_bound, std::min, std::max
#include <cstring>    // for std::memcpy

// Loads in and models an exact tabulated shape function, interpolating for points in between.

extern "C" {
    OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model);
    OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name);
}

class TabulatedModel : public OblModelBase {
public:
    TabulatedModel(double req_nounits, double mass_nounits, double omega_nounits, int model);
    TabulatedModel(double req_nounits, double mass_nounits, double omega_nounits, std::string model_name);
    double R_at_costheta(const double& costheta) const override;
    double Dtheta_R(const double& costheta);

private:
    std::map<double, double> mu_radii_map; // Store mu values and corresponding radii in nounits
    void loadShapeFunction(int model);
    void loadShapeFunction(std::string model);
};
#endif // POLYOBLMODELBASE_H
