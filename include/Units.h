// Units.h
//
// This files contains the definitions of quantities used for
// converting dimensionful quantities to dimensionless quantities
// and declares some functions useful for doing so.
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


#ifndef UNITS_H
#define UNITS_H

// Some useful values a,b,c:
// Mass density = energy density/c^2 = {-3, 1, 0}
// Pressure = {-1,1,-2}
// Enthalpy = {2,0,-2}
// Number density = {-3,0,0}

namespace Units {
  const double GMC2 = 1.47662512; // G x MSUN/C^2 in km (2017)
  const double C = 2.99792458e10; // speed of light in vacuo (NIST recommended, 19 Oct 2004, 2002 CODATA)
  const double G = 6.6742e-8;  // gravitational constant (NIST recommended, 19 Oct 2004, 2002 CODATA) 
  const double KAPPA = (1.0e-15*C*C/G); // length scale  (consider this as some number expressed in cm^2)
  const double MSUN = 1.9891e33; // ("IAU Style Manual", Wilkinson, G.A. Comm. 5, IAU Transactions XXB (1987)) 
                                        // solar mass
  const double MBARYON = 1.67e-24;  // baryon mass  
  const double PI = 3.14159265358979323846; // (Eric W. Weisstein. "Pi." From MathWorld--A Wolfram Web Resource. 
                                             //    http://mathworld.wolfram.com/Pi.html)
  const double H_PLANCK = 6.6260693e-27; // erg-seconds (CODATA 2002).  This isn't hbar.  hbar = H_PLANCK/(2*PI)
  const double EV = 1.60217653e-12; // 1 eV in ergs (CODATA 2002)
  const double K_BOLTZ = 1.3806505e-16; // erg / Kelvin  (CODATA 2002)
  const double PARSEC = 3.08572e18; // in cm (agrees with IAU, above)
  const double STEPH_BOLTZMANN = 5.670374419e-5; //In ergs/cm^2/k^4, Codata (2022) https://physics.nist.gov/cuu/Constants/Table/allascii.txt

  const double EV_TO_K = 1.160451812e4; // Codata (2022) https://www.physics.nist.gov/cgi-bin/cuu/Value?evk, conversion factor between ev and kelvin.
  const double K_TO_EV = 8.617333262e-5; // Codata (2022) https://physics.nist.gov/cuu/Constants/Table/allascii.txt 
  const double EV_TO_HZ = 2.417989242e14; // Codata (2022) https://physics.nist.gov/cuu/Constants/Table/allascii.txt 

  struct Dimensions{ int a, b, c; };
  const Dimensions MASS = {0, 1, 0};
  const Dimensions LENGTH = {1, 0, 0};
  const Dimensions TIME = {0, 0, 1};
  const Dimensions MASSPERLENGTH = {-1, 1, 0}; 
  const Dimensions INVLENGTH = {-1, 0, 0}; 
  const Dimensions INVTIME = {0, 0, -1};
  const Dimensions MASSDENSITY = {-3, 1, 0};
  const Dimensions PRESSURE = {-1,1,-2};
  const Dimensions ENTHALPY = {2,0,-2};
  const Dimensions NUMBERDENSITY = {-3,0,0};
  const Dimensions MOMINERTIA = {2, 1, 0};
  const Dimensions ANGMOMENTUM = {2, 1, -1};
  const Dimensions ACCEL = {1, 0, -2};
  const Dimensions VELOCITY = {1, 0, -1};
  const Dimensions ENERGY = {2, 1, -2};
  const Dimensions ENERGYFLUX = {0, 1, -3}; // erg s^-1 cm^-2
  const Dimensions INTENSITY = {0, 1, -2}; // erg s^-1 cm^-2 Hz^-1
  double cgs_to_nounits(const double& cgs, const Dimensions& d);
  double nounits_to_cgs(const double& nounits, const Dimensions& d);
}

#endif // UNITS_H

