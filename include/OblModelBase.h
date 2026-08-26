// OblModelBase.h
//
// Interface describing an Oblateness model
// (C) Coire Cadeau, 2007

// Source (C) Coire Cadeau 2007, all rights reserved.
//
// Permission is granted for private use only, and not
// distribution, either verbatim or of derivative works,
// in whole or in part.
//
// The code is not thoroughly tested or guaranteed for
// any particular use.

#ifndef OBLMODELBASE_H
#define OBLMODELBASE_H

#include <exception>
#include <string>

class OblModelBase {
 protected:
  std::string model_name;
  double req_nounits, mass_nounits, omega_nounits;
  OblModelBase(double req_nounits, double mass_nounits, double omega_nounits) : req_nounits(req_nounits), mass_nounits(mass_nounits), omega_nounits(omega_nounits) {}
  OblModelBase(double req_nounits, double mass_nounits, double omega_nounits, std::string model_name) : req_nounits(req_nounits), mass_nounits(mass_nounits), omega_nounits(omega_nounits), model_name(model_name) {}

 public:
  virtual ~OblModelBase() { }
  virtual double R_at_costheta( const double& costheta ) const = 0;
  virtual double Dtheta_R( const double& costheta ) = 0;

  double Req_nounits() const;
  double Mass_nounits() const;
  double Omega_nounits() const;
  void get_model_name(std::string* model) const;

  double cos_gamma(const double& costheta) ; // return cosine of angle between radial and normal

  private:
   double z(const double& costheta) const;
   double f(const double& costheta) ;
   double get_zeta() const;
};

#endif // OBLMODELBASE_H
