#ifndef BAUBOCKETALHT_H
#define BAUBOCKETALHT_H

#include "../../include/OblModelBase.h"
#include <exception>
#include <string>

//Encodes the Baubock Et Al 2013 Paper's Shape Function
// RELATIONS BETWEEN NEUTRON-STAR PARAMETERS IN THE HARTLE–THORNE APPROXIMATION
// This one remains within HT, does not go to BL coordinates.
//See: https://iopscience.iop.org/article/10.1088/0004-637X/777/1/68/pdf

extern "C" {
    OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model);
    OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name);
}

class BaubockOblModel : public OblModelBase {
 protected:

 public:
  BaubockOblModel(double req_nounits, double mass_nounits, double omega_nounits);

  double R_at_costheta( const double& costheta ) const ;
  double Dtheta_R( const double& costheta ) ;

 private:  
  double Req_HT;
  double findReqHT(double req_nounits);

};
#endif 
