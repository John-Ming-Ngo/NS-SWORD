#ifndef SPYY_PAPIGKIOTIS_RP_H
#define SPYY_PAPIGKIOTIS_RP_H

#include "../../include/OblModelBase.h"
#include <string>

extern "C" {
OblModelBase* createOblModelBase(const double&, const double&, const double&, int);
OblModelBase* createOblModelBaseStr(const double&, const double&, const double&, std::string);
}

class SPYYPapigkiotisRp : public OblModelBase {
 public:
  SPYYPapigkiotisRp(double req, double mass, double omega, int model = 0);
  SPYYPapigkiotisRp(double req, double mass, double omega, std::string model_name);

  double R_at_costheta(const double& costheta) const;
  double Dtheta_R(const double& costheta);

  static double compactness(double mass, double req);
  static double spin(double omega, double mass, double req);
  double polar_ratio() const;

 private:
  void load_spyy_coefficients();
  double evaluate_coefficient(const double c[7]) const;
  double g(double mu) const;
  double dg_dmu(double mu) const;
  double validated_radical(double mu) const;

  double C_, sigma_, polar_ratio_, e2_, a2_, a4_;
  double a2_coeff_[7], a4_coeff_[7];
};

#endif
