#ifndef PAPIGKIOTIS_H
#define PAPIGKIOTIS_H

#include "../../include/OblModelBase.h"
#include <map>    
#include <cmath>       

#include <exception>
#include <string>

extern "C" {
    OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model);
    OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name);
}

class Papigkiotis : public OblModelBase {
 protected:

 public:
  std::map<double, double> mu_radii_map;
  std::map<double, double> mu_radii_derivative_map;

  Papigkiotis(double req_nounits, double mass_nounits, double omega_nounits, int model);
  Papigkiotis(double req_nounits, double mass_nounits, double omega_nounits, std::string model_name);

  static double C( const double& Mass_nounits, const double& Req_nounits );
  static double sigma( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits );

  double R_at_costheta( const double& costheta ) const ;
  double Dtheta_R( const double& costheta ) ;
  
 private:  

};
#endif // PAPIGKIOTIS_H
