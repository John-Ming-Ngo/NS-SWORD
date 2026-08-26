// Morsink.h
//
// Interface describing a Polynomial Oblateness model
// (C) Coire Cadeau, 2007

// Source (C) Coire Cadeau 2007, all rights reserved.
//
// Permission is granted for private use only, and not
// distribution, either verbatim or of derivative works,
// in whole or in part.
//
// The code is not thoroughly tested or guaranteed for
// any particular use.

// SMM: Nov 30, 2009
// Added Rspot_nounits as a parameter and function.

// Encodes the Morsink et al 2007 Shape Function https://iopscience.iop.org/article/10.1086/518648/meta
// THE OBLATE SCHWARZSCHILD APPROXIMATION FOR LIGHT CURVES OF RAPIDLY ROTATING NEUTRON STARS

#ifndef MORSINK_H
#define MORSINK_H

#include "../../include/OblModelBase.h"
#include <exception>
#include <string>

extern "C" {
    OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model);
    OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name);
}

class Morsink : public OblModelBase {
 private:
  int model;
  double a_0_0, a_0_1, a_0_2, a_2_0, a_2_1, a_2_2, a_4_0, a_4_1, a_4_2;
  double zeta, eps;

 public:
  Morsink(const double& req_nounits, const double& mass_nounits, const double& omega_nounits);
  Morsink(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model);
  Morsink(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name);

  double R_at_costheta( const double& costheta ) const ;
  double Dtheta_R( const double& costheta ) ;

  static double zetaparam( const double& Mass_nounits, const double& Req_nounits );
  static double epsparam( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits );
  virtual ~Morsink() { }
  static bool validModel(int model);

  double a0() const;
  double a2() const;
  double a4() const;

 protected:  
  double get_zeta() const;
  double get_eps() const;
};

#endif // MORSINK_H
