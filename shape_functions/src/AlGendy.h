#ifndef ALGENDY_H
#define ALGENDY_H

#include "../../include/OblModelBase.h"
#include <exception>
#include <string>

// Encodes the 2014 AlGendy and Morsink 2014 paper's shape function which is used with NICER's analysis
//  UNIVERSALITY OF THE ACCELERATION DUE TO GRAVITY ON THE SURFACE OF A RAPIDLY ROTATING NEUTRON STAR
// See: https://iopscience.iop.org/article/10.1088/0004-637X/791/2/78/pdf

extern "C" {
    OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model);
    OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name);
}

class AlGendyOblModel : public OblModelBase {
 public:
  AlGendyOblModel(double req_nounits, double mass_nounits, double omega_nounits, int model);
  AlGendyOblModel(double req_nounits, double mass_nounits, double omega_nounits, std::string model_name_input);

  double R_at_costheta( const double& costheta ) const ;
  double Dtheta_R( const double& costheta ) ;

  static double zetaparam( const double& Mass_nounits, const double& Req_nounits );
  static double epsparam( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits );
  static bool validModel(int model);

  double a2() const;
  
 private:
   int model;
   double zeta, eps;

   double c_2_0, c_2_1 = 0;

 protected:  
  double get_zeta() const;
  double get_eps() const;
};
#endif // ALGENDY_H
