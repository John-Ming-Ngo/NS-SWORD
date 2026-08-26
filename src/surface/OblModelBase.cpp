#include "../../include/OblModelBase.h"
#include <cmath>

double OblModelBase::Req_nounits() const {
    return req_nounits;
}
double OblModelBase::Mass_nounits() const {
    return mass_nounits;
}
double OblModelBase::Omega_nounits() const {
    return omega_nounits;
}
void OblModelBase::get_model_name(std::string* model) const {
    *model = model_name;
}

double OblModelBase::get_zeta() const {
    return double( Mass_nounits() / Req_nounits() );
}

double OblModelBase::z(const double& costheta) const {
    return 1.0 / sqrt(1.0 - 2.0 * get_zeta() * Req_nounits() / R_at_costheta(costheta)) - 1.0;
}

double OblModelBase::f(const double& costheta)  {
    return (1.0 + z(costheta)) * Dtheta_R(costheta) / R_at_costheta(costheta);
}

double OblModelBase::cos_gamma(const double& costheta)  {
    return 1.0 / sqrt(1.0 + pow(f(costheta), 2.0));
}