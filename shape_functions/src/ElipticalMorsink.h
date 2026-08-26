// ElipticalMorsink.h
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

#ifndef ELIPTICALMORSINK_H
#define ELIPTICALMORSINK_H

#include "../../include/OblModelBase.h"
#include <exception>
#include <string>

// Encodes the Baubock Et Al. 2012 Shape function which utilizes a modified enforced eliptical version of the Morsink et al 2007 shape.
// See: https://iopscience.iop.org/article/10.1088/0004-637X/753/2/175/pdf
//  A RAY-TRACING ALGORITHM FOR SPINNING COMPACT OBJECT SPACETIMES WITH ARBITRARY QUADRUPOLE MOMENTS. II.NEUTRONSTARS

extern "C" {
    OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model);
    OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name);
}

class ElipticalMorsink : public OblModelBase {
 private:
  int model;
  double a_0_0, a_0_1, a_0_2, a_2_0, a_2_1, a_2_2, a_4_0, a_4_1, a_4_2;
  double zeta, eps;

 public:
  ElipticalMorsink(const double& req_nounits, const double& mass_nounits, const double& omega_nounits);
  ElipticalMorsink(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model);
  ElipticalMorsink(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name);

  double R_at_costheta( const double& costheta ) const ;
  double Dtheta_R( const double& costheta ) ;

  static double zetaparam( const double& Mass_nounits, const double& Req_nounits );
  static double epsparam( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits );
  virtual ~ElipticalMorsink() { }
  static bool validModel(int model);

  double a0() const;
  double a2() const;
  double a4() const;

 protected:  
  double get_zeta() const;
  double get_eps() const;
};

#endif // ELIPTICALMORSINK_H
