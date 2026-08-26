#ifndef SILVAETAL_H
#define SILVAETAL_H

#include "../../include/OblModelBase.h"
#include <exception>
#include <string>

// Models the eliptical shape function given by Silva et al 2021 in: https://ui.adsabs.harvard.edu/abs/2021PhRvD.103f3038S/abstract
// Surface of rapidly-rotating neutron stars: Implications to neutron star parameter estimation

extern "C" {
    OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model);
    OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name);
}

class SilvaOblModel : public OblModelBase {
 public:
  SilvaOblModel(double req_nounits, double mass_nounits, double omega_nounits);
  SilvaOblModel(double req_nounits, double mass_nounits, double omega_nounits, int model);
  SilvaOblModel(double req_nounits, double mass_nounits, double omega_nounits, std::string model_name);

  double R_at_costheta( const double& costheta ) const ;
  double Dtheta_R( const double& costheta ) ;

  static double zetaparam( const double& Mass_nounits, const double& Req_nounits );
  static double epsparam( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits );
  static bool validModel(int model);

  double y(double c_0_0, double c_1_2_0, double c_1_0, double c_0_1, double c_1_1, double c_2_0, double c_0_2, double eps, double zeta) const;

  double eccentricity() const;
  double a2() const;
  double a4() const;

  double g(const double& costheta) const;
  double dg_dcostheta(const double& costheta) const;
  
 private:
   int model;
   double zeta, eps;

   double e_c_0_0, e_c_1_2_0, e_c_1_0, e_c_0_1, e_c_1_1, e_c_2_0, e_c_0_2, 
   a2_c_0_0, a2_c_1_2_0, a2_c_1_0, a2_c_0_1, a2_c_1_1, a2_c_2_0, a2_c_0_2, 
   a4_c_0_0, a4_c_1_2_0, a4_c_1_0, a4_c_0_1, a4_c_1_1, a4_c_2_0, a4_c_0_2 = 0;

 protected:  
  double get_zeta() const;
  double get_eps() const;
};
#endif // POLYOBLMODELBASE_H
