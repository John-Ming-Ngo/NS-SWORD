// Units.cpp
//
// This file contains the source code for functions that convert
// cgs units to dimensionless units used by this code, and vice-versa.
//
// If you need units of L^a M^b T^c, you have to multiply by
// KAPPA^x G^y C^z, where:
//
// x = (a+b+c)/2
// y = 2b - c
// z = -b
//
// and divide the dimensionful quantities by the same factor to get the
// dimensionless version.
// These functions actually take an instance of the Dimensions struct defined
// in the header (which has int members a, b, c)
//
// (C) Coire Cadeau, 2007

// Source (C) Coire Cadeau 2007, all rights reserved.
//
// Permission is granted for private use only, and not
// distribution, either verbatim or of derivative works,
// in whole or in part.
//
// The code is not thoroughly tested or guaranteed for
// any particular use.

#include "../../include/Units.h"
#include <cmath>
// #include <iostream>

using namespace Units;

double Units::cgs_to_nounits(const double& cgs, const Dimensions& d) {
  // the ".0" on the numbers below are crucial.  thank God for long weekends.
  double factor( pow(KAPPA,(d.a+d.b+d.c)/2.0)*pow(G,-d.b)*pow(C,2.0*d.b-d.c) );

  /*
  if(d.a == Units::LENGTH.a && d.b == Units::LENGTH.b && d.c == Units::LENGTH.c) {
    std::cerr << "Length scale conversion cgs->nounits requested: " << std::endl
	      << "Dimensions a b c" << d.a << " " << d.b << " " << d.c << std::endl
	      << "exponent: " << (d.a+d.b+d.c)/2 << std::endl
	      << "sqrt(KAPPA) = " << sqrt(KAPPA) << std::endl
	      << "sqrt(1.0e-15*C*C/G) = " << sqrt(1.0e-15*C*C/G) << std::endl
	      << "factor = " << factor << std::endl
	      << "cgs = " << cgs << std::endl
	      << "cgs/factor = " << cgs/factor << std::endl;
	      }*/

  return cgs/factor;
}

double Units::nounits_to_cgs(const double& nounits, const Dimensions& d){
  double factor( pow(KAPPA,(d.a+d.b+d.c)/2.0)*pow(G,-d.b)*pow(C,2.0*d.b-d.c) );
  return nounits*factor;
}

