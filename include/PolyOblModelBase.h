// PolyOblModelBase.h
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

#ifndef POLYOBLMODELBASE_H
#define POLYOBLMODELBASE_H

#include "OblModelBase.h"
#include <exception>
#include <string>

class PolyOblModelBase : public OblModelBase {
 private:
  int model;
  double a_0_0, a_0_1, a_0_2, a_2_0, a_2_1, a_2_2, a_4_0, a_4_1, a_4_2;
  double zeta, eps;

 public:
  PolyOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits);
  PolyOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model);

  double R_at_costheta( const double& costheta ) ;
  double Dtheta_R( const double& costheta ) ;

  static double zetaparam( const double& Mass_nounits, const double& Req_nounits );
  static double epsparam( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits );
  virtual ~PolyOblModelBase() { }
  static bool validModel(int model);

  double a0() const;
  double a2() const;
  double a4() const;

 protected:  
  double get_zeta() const;
  double get_eps() const;
};

#endif // POLYOBLMODELBASE_H
