#define _USE_MATH_DEFINES

#include "BaubockEtAlBL.h"

#include "../../include/Units.h"
#include <iostream>
#include <cstdlib>  // For std::atof
#include <cmath>
#include <stdexcept>
#include <math.h>


BaubockOblModel::BaubockOblModel(double req_nounits, double mass_nounits, double omega_nounits)
: OblModelBase(req_nounits, mass_nounits, omega_nounits), Req_HT(req_nounits), radius_scale(1.0) {
  const double raw_equator = R_at_costheta(0.0);
  const double candidate_scale = req_nounits / raw_equator;
  bool valid_profile = std::isfinite(candidate_scale) && candidate_scale > 0.0;
  const double dtheta = (Units::PI - 1.0e-6) / 128.0;
  for (int i = 0; valid_profile && i < 128; i++) {
    const double theta = i * dtheta + 1.0e-6 + 0.5 * dtheta;
    const double mu = cos(theta);
    const double radius = R_at_costheta(mu) * candidate_scale;
    valid_profile = std::isfinite(radius) && radius > 0.0;
  }
  if (valid_profile) {
    radius_scale = candidate_scale;
  }
  else {
    Req_HT = findReqHT(req_nounits);
  }
  model_name = "Baubock_BL";
}

// Constants
const double G = Units::G * 1e-3;  // Gravitational constant in m^3/kg/s^2
const double c = Units::C / 100;        // Speed of light in m/s
double BaubockOblModel::findReqHT(double req_nounits) {
  Req_HT = req_nounits;
  double factor = 0.2;
  const double epsilon = 1.0e-8;
  int iterations = 0;

  while (req_nounits - R_at_costheta(0.0) > epsilon) {
    Req_HT *= 1.0 + factor;
    if (req_nounits - R_at_costheta(0.0) < 0.0) {
      Req_HT /= 1.0 + factor;
      factor /= 2.0;
    }
    if (!std::isfinite(Req_HT) || ++iterations > 10000) {
      throw std::runtime_error("Baubock: equatorial-radius search did not converge");
    }
  }
  while (R_at_costheta(0.0) - req_nounits > epsilon) {
    Req_HT *= 1.0 - factor;
    if (R_at_costheta(0.0) - req_nounits < 0.0) {
      Req_HT /= 1.0 - factor;
      factor /= 2.0;
    }
    if (!std::isfinite(Req_HT) || ++iterations > 10000) {
      throw std::runtime_error("Baubock: equatorial-radius search did not converge");
    }
  }
  return Req_HT;
}

double zeta(double M_0, double R_0) {
  return (G * M_0) / (R_0 * pow(c, 2));
}

// Equation (1): Second-order Legendre Polynomial P_2
double P_2(double cos_theta) {
  return 0.5 * (3 * pow(cos_theta, 2) - 1);
}

// Equation (2): epsilon_0 definition
double epsilon_0(double f, double f_0) {
  return f / f_0;
}

// Equation (3): f_0 calculation
double f_0(double M_0, double R_0) {
  return sqrt(G * M_0 / pow(R_0, 3)) / (2*M_PI); //The original equation does not contain the division by 2pi, but this appears to be a dimensionless to dimensionful error.
}

// Equation (4): epsilon_s using epsilon_0 and epsilon_s_star
double epsilon_s(double epsilon_0, double epsilon_s_star) {
  return pow(epsilon_0, 2) * epsilon_s_star;
}

// Equation (5): a_star calculation
double a_star(double zeta) {
  return 1.1035 - 2.146 * zeta + 4.5756 * pow(zeta, 2);
}

// Equation (6): Q_bar calculation
double Q_bar(double q_star, double a_star) {
  return q_star / pow(a_star, 2);
}

// Equation (7): I_bar calculation
double I_bar(double a_star, double zeta) {
  return a_star * pow(zeta, -1.5);
}

// Equation (8): ln(Q_bar) calculation
double ln_Q_bar(double ln_I_bar) {
  return -2.014 + 0.601 * ln_I_bar + 1.10 * pow(ln_I_bar, 2) - 0.412 * pow(ln_I_bar, 3) + 0.0459 * pow(ln_I_bar, 4);
}

// Equation (9): epsilon_s final form involving xi_2
double epsilon_s_final(double R, double xi_2) {
  return -(3.0 / (2.0 * R)) * xi_2;
}

// Equation (10): Complex epsilon_s_star calculation
double epsilon_s_star(double zeta, double a_star, double q_star) {
  if (zeta == 0.5) {
    throw std::invalid_argument("Logarithmic singularity at zeta = 0.5");
  }

  double part1 = 8 * pow(zeta, 2) - 32 * a_star * pow(zeta, 7.0 / 2.0);
  double part2 = (pow(a_star, 2) - q_star) * (45 - 135 * zeta + 60 * pow(zeta, 2) + 30 * pow(zeta, 3));
  double part3 = 24 * pow(a_star, 2) * pow(zeta, 4) + 8 * pow(a_star, 2) * pow(zeta, 5) - 48 * pow(a_star, 2) * pow(zeta, 6);
  double part4 = 45 * (pow(a_star, 2) - q_star) * pow(1 - 2 * zeta, 2) * log(1 - 2 * zeta);
    
  return (1.0 / (32 * pow(zeta, 3))) * (2 * zeta * (part1 + part2 + part3) + part4);
}

// Equation (11): R_BL calculation
double R_BL(double R_HT, double theta, double M, double a) {
  double GM_c2 = G * M / pow(c, 2);
  double R_HT_cubed = pow(R_HT, 3);
   
  double factor = pow(GM_c2 * a, 2) / (2 * R_HT_cubed);
  double part1 = (R_HT + 2 * GM_c2) * (R_HT - GM_c2);
  double part2 = pow(cos(theta), 2) * (R_HT - 2 * GM_c2) * (R_HT + 3 * GM_c2);
    
  return R_HT - factor * (part1 - part2);
}

// Equation (12): Final equation for R(theta)
double R_theta(double R_0, double xi_2, double theta) {
    return R_0 + xi_2 * P_2(cos(theta));
}

double dR_dtheta(double R_0, double xi_2, double theta) {
    return -3*xi_2*cos(theta)*sin(theta);
}

double R_0_BL(double M_0, double R_0, double f, double theta) {
  double zeta_val = zeta(M_0, R_0);
  //std::cout << "zeta(M_0, R_0): " << zeta_val << std::endl;
  double f_0_val = f_0(M_0, R_0);
  //std::cout << "f_0(M_0, R_0): " << f_0_val << std::endl;
  double epsilon_0_val = epsilon_0(f, f_0_val);
  //std::cout << "epsilon_0(f, f_0_val): " << epsilon_0_val << std::endl;

  // Calculate intermediate values
  double a_star_val = a_star(zeta_val);
  //std::cout << "a_star(zeta_val): " << a_star_val << std::endl;
  double a_val = epsilon_0_val * a_star_val;
  //std::cout << "a(a_star, epsilon_0): " << a_val << std::endl;
  double I_bar_val = I_bar(a_star_val, zeta_val);
  //std::cout << "I_bar(a_star_val, zeta_val): " << I_bar_val << std::endl;
  double lnQ_val = ln_Q_bar(log(I_bar_val));
  //std::cout << "ln_Q_bar(log(I_bar_val)): " << lnQ_val << std::endl;
  double Q_bar_val = exp(lnQ_val);
  //std::cout << "Q_bar: " << Q_bar_val << std::endl;
  double q_star_val = Q_bar_val * pow(a_star_val, 2);
  //std::cout << "q_star(Q_bar_val, a_star_val): " << q_star_val << std::endl;
  double epsilon_s_star_val = epsilon_s_star(zeta_val, a_star_val, q_star_val);
  //std::cout << "epsilon_s_star(zeta_val, a_star_val, q_star_val): " << epsilon_s_star_val << std::endl;
  double epsilon_s_val = epsilon_s(epsilon_0_val, epsilon_s_star_val);
  //std::cout << "epsilon_s(epsilon_0_val, epsilon_s_star_val): " << epsilon_s_val << std::endl;
  double xi_2_val = -(2.0 * R_0 / 3.0) * epsilon_s_val;
  //std::cout << "xi_2_val: " << xi_2_val << std::endl;

  // Calculate R(theta)
  double R_theta_val = R_theta(R_0, xi_2_val, theta);
  //std::cout << "R_theta(R_0, xi_2_val, theta): " << R_theta_val << std::endl;
  double R_theta_BL_val = R_BL(R_0, theta, M_0, a_val);
  //std::cout << "R_BL(R_theta_val, theta, M_0): " << R_theta_BL_val << std::endl;

  return R_theta_BL_val;
}

double R_theta_BL(double M_0, double R_0, double f, double theta) {
  double zeta_val = zeta(M_0, R_0);
  //std::cout << "zeta(M_0, R_0): " << zeta_val << std::endl;
  double f_0_val = f_0(M_0, R_0);
  //std::cout << "f_0(M_0, R_0): " << f_0_val << std::endl;
  double epsilon_0_val = epsilon_0(f, f_0_val);
  //std::cout << "epsilon_0(f, f_0_val): " << epsilon_0_val << std::endl;
  //std::cout << "GM/R^3: " << G*M_0/pow(R_0, 3) << std::endl;

  // Calculate intermediate values
  double a_star_val = a_star(zeta_val);
  //std::cout << "a_star(zeta_val): " << a_star_val << std::endl;
  double a_val = epsilon_0_val * a_star_val;
  //std::cout << "a(a_star, epsilon_0): " << a_val << std::endl;
  double I_bar_val = I_bar(a_star_val, zeta_val);
  //std::cout << "I_bar(a_star_val, zeta_val): " << I_bar_val << std::endl;
  double lnQ_val = ln_Q_bar(log(I_bar_val));
  //std::cout << "ln_Q_bar(log(I_bar_val)): " << lnQ_val << std::endl;
  double Q_bar_val = exp(lnQ_val);
  //std::cout << "Q_bar: " << Q_bar_val << std::endl;
  double q_star_val = Q_bar_val * pow(a_star_val, 2);
  //std::cout << "q_star(Q_bar_val, a_star_val): " << q_star_val << std::endl;
  double epsilon_s_star_val = epsilon_s_star(zeta_val, a_star_val, q_star_val);
  //std::cout << "epsilon_s_star(zeta_val, a_star_val, q_star_val): " << epsilon_s_star_val << std::endl;
  double epsilon_s_val = epsilon_s(epsilon_0_val, epsilon_s_star_val);
  //std::cout << "epsilon_s(epsilon_0_val, epsilon_s_star_val): " << epsilon_s_val << std::endl;
  double xi_2_val = -(2.0 * R_0 / 3.0) * epsilon_s_val;
  //std::cout << "xi_2_val: " << xi_2_val << std::endl;

  // Calculate R(theta)
  double R_theta_val = R_theta(R_0, xi_2_val, theta);
  //std::cout << "R_theta(R_0, xi_2_val, theta): " << R_theta_val << std::endl;
  double R_theta_BL_val = R_BL(R_theta_val, theta, M_0, a_val);
  //std::cout << "R_BL(R_theta_val, theta, M_0): " << R_theta_BL_val << std::endl;

  return R_theta_BL_val;
}

double R_theta_HT(double M_0, double R_0, double f, double theta) {
  double zeta_val = zeta(M_0, R_0);
  //std::cout << "zeta(M_0, R_0): " << zeta_val << std::endl;
  double f_0_val = f_0(M_0, R_0);
  //std::cout << "f_0(M_0, R_0): " << f_0_val << std::endl;
  double epsilon_0_val = epsilon_0(f, f_0_val);
  //std::cout << "epsilon_0(f, f_0_val): " << epsilon_0_val << std::endl;

  //std::cout << "GM/R^3: " << G*M_0/pow(R_0, 3) << std::endl;

  // Calculate intermediate values
  double a_star_val = a_star(zeta_val);
  //std::cout << "a_star(zeta_val): " << a_star_val << std::endl;
  double a_val = epsilon_0_val * a_star_val;
  //std::cout << "a(a_star, epsilon_0): " << a_val << std::endl;
  double I_bar_val = I_bar(a_star_val, zeta_val);
  //std::cout << "I_bar(a_star_val, zeta_val): " << I_bar_val << std::endl;
  double lnQ_val = ln_Q_bar(log(I_bar_val));
  //std::cout << "ln_Q_bar(log(I_bar_val)): " << lnQ_val << std::endl;
  double Q_bar_val = exp(lnQ_val);
  //std::cout << "Q_bar: " << Q_bar_val << std::endl;
  double q_star_val = Q_bar_val * pow(a_star_val, 2);
  //std::cout << "q_star(Q_bar_val, a_star_val): " << q_star_val << std::endl;
  double epsilon_s_star_val = epsilon_s_star(zeta_val, a_star_val, q_star_val);
  //std::cout << "epsilon_s_star(zeta_val, a_star_val, q_star_val): " << epsilon_s_star_val << std::endl;
  double epsilon_s_val = epsilon_s(epsilon_0_val, epsilon_s_star_val);
  //std::cout << "epsilon_s(epsilon_0_val, epsilon_s_star_val): " << epsilon_s_val << std::endl;
  double xi_2_val = -(2.0 * R_0 / 3.0) * epsilon_s_val;
  //std::cout << "xi_2_val: " << xi_2_val << std::endl;

  // Calculate R(theta)
  double R_theta_val = R_theta(R_0, xi_2_val, theta);
  //std::cout << "R_theta(R_0, xi_2_val, theta): " << R_theta_val << std::endl;

  return R_theta_val;
}

double dR_theta_HT(double M_0, double R_0, double f, double theta) {
  double zeta_val = zeta(M_0, R_0);
  //std::cout << "zeta(M_0, R_0): " << zeta_val << std::endl;
  double f_0_val = f_0(M_0, R_0);
  //std::cout << "f_0(M_0, R_0): " << f_0_val << std::endl;
  double epsilon_0_val = epsilon_0(f, f_0_val);
  //std::cout << "epsilon_0(f, f_0_val): " << epsilon_0_val << std::endl;

  //std::cout << "GM/R^3: " << G*M_0/pow(R_0, 3) << std::endl;

  // Calculate intermediate values
  double a_star_val = a_star(zeta_val);
  //std::cout << "a_star(zeta_val): " << a_star_val << std::endl;
  double a_val = epsilon_0_val * a_star_val;
  //std::cout << "a(a_star, epsilon_0): " << a_val << std::endl;
  double I_bar_val = I_bar(a_star_val, zeta_val);
  //std::cout << "I_bar(a_star_val, zeta_val): " << I_bar_val << std::endl;
  double lnQ_val = ln_Q_bar(log(I_bar_val));
  //std::cout << "ln_Q_bar(log(I_bar_val)): " << lnQ_val << std::endl;
  double Q_bar_val = exp(lnQ_val);
  //std::cout << "Q_bar: " << Q_bar_val << std::endl;
  double q_star_val = Q_bar_val * pow(a_star_val, 2);
  //std::cout << "q_star(Q_bar_val, a_star_val): " << q_star_val << std::endl;
  double epsilon_s_star_val = epsilon_s_star(zeta_val, a_star_val, q_star_val);
  //std::cout << "epsilon_s_star(zeta_val, a_star_val, q_star_val): " << epsilon_s_star_val << std::endl;
  double epsilon_s_val = epsilon_s(epsilon_0_val, epsilon_s_star_val);
  //std::cout << "epsilon_s(epsilon_0_val, epsilon_s_star_val): " << epsilon_s_val << std::endl;
  double xi_2_val = -(2.0 * R_0 / 3.0) * epsilon_s_val;
  //std::cout << "xi_2_val: " << xi_2_val << std::endl;

  // Calculate R(theta)
  double dR_theta_val = dR_dtheta(R_0, xi_2_val, theta);
  //std::cout << "R_theta(R_0, xi_2_val, theta): " << R_theta_val << std::endl;

  return dR_theta_val;
}

double internal_function(double M_0, double R_0, double f, double theta) {
  double radius_HT = R_theta_HT(M_0, R_0, f, theta);
  double GM_c2 = G * M_0 / pow(c, 2);
  return ((radius_HT + 2*GM_c2) * (radius_HT-GM_c2)) - ((pow(cos(theta), 2) * (radius_HT - 2 * GM_c2) * (radius_HT + 3 *GM_c2)));
}

double d_internal_function(double M_0, double R_0, double f, double theta) {
  double dradius_HT = dR_theta_HT(M_0, R_0, f, theta);
  double radius_HT = R_theta_HT(M_0, R_0, f, theta);
  double GM_c2 = G * M_0 / pow(c, 2);

  return (dradius_HT*(2*radius_HT+GM_c2))
  - 2*cos(theta)*(-1*sin(theta))*(radius_HT-2*GM_c2)*(radius_HT+3*GM_c2)
  - pow(cos(theta), 2) * dradius_HT * (2*radius_HT+GM_c2);
}

double dR_dtheta_BL(double M_0, double R_0, double f, double theta) {
  double zeta_val = zeta(M_0, R_0);
  double f_0_val = f_0(M_0, R_0);
  double epsilon_0_val = epsilon_0(f, f_0_val);

  // Calculate intermediate values
  double a_star_val = a_star(zeta_val);
  double a_val = epsilon_0_val * a_star_val;
  double I_bar_val = I_bar(a_star_val, zeta_val);
  double lnQ_val = ln_Q_bar(log(I_bar_val));
  double Q_bar_val = exp(lnQ_val);
  double q_star_val = Q_bar_val * pow(a_star_val, 2);
  double epsilon_s_star_val = epsilon_s_star(zeta_val, a_star_val, q_star_val);
  double epsilon_s_val = epsilon_s(epsilon_0_val, epsilon_s_star_val);
  double xi_2_val = -(2.0 * R_0/3.0) * epsilon_s_val;

  // Calculate R(theta)
  double R_theta_val = R_theta(R_0, xi_2_val, theta);
  double R_theta_BL_val = R_BL(R_theta_val, theta, M_0, a_val);

  // Other important values
  double dradius_HT = dR_theta_HT(M_0, R_0, f, theta);
  double radius_HT = R_theta_HT(M_0, R_0, f, theta);
  double GM_c2 = G * M_0 / pow(c, 2);
  double internal_function_val = internal_function(M_0, R_0, f, theta);
  double d_internal_function_val = d_internal_function(M_0, R_0, f, theta);

  //Calculate the derivative

  return dradius_HT
  - (pow(a_val * GM_c2, 2))/2 * ((-3/pow(radius_HT, 4) * dradius_HT * internal_function_val + 1/pow(radius_HT, 3) * d_internal_function_val));
}

double BaubockOblModel::R_at_costheta( const double& costheta ) const {
  // Return R(theta) in "nounits".
  // note that the user supplies cos(theta) and not theta.

  //Step 0: Get m, f, req, and theta.
  double mass_unitless = Mass_nounits();
  double omega_unitless = Omega_nounits();
  double radius_unitless = Req_HT;
  double theta = acos(costheta);

  //Step 1: Turn unitless to cgs.
  double mass_cgs = Units::nounits_to_cgs(mass_unitless, Units::MASS);
  double f_cgs = Units::nounits_to_cgs(omega_unitless, Units::INVTIME)/(2.0*Units::PI);
  double radius_cgs = Units::nounits_to_cgs(radius_unitless, Units::LENGTH);

  //Step 2: Turn cgs to mks.
  double mass_mks = mass_cgs/1000;
  double f_mks = f_cgs;
  double radius_mks = radius_cgs/100;

  //Step 3: Use R_theta_BL(double M_0, double R_0, double f, double theta) to get the radii.
  double radius_theta_bl = R_theta_BL(mass_mks, radius_mks, f_mks, theta);
  double radius_equator_bl = R_theta_BL(mass_mks, radius_mks, f_mks, Units::PI/2);

  //Step 4: Turn the value to cgs, and then to unitless.
  double radius_theta_bl_cgs = radius_theta_bl * 100;
  double radius_theta_bl_unitless = Units::cgs_to_nounits(radius_theta_bl_cgs, Units::LENGTH);

  //Error diagnostic
  //std::cout 
  //  << "costheta: " << costheta 
  //  << ", theta: " << acos(costheta) 
  //  << ", radius_mks: " << radius_mks 
  //  << ", radius calculated: " << radius_theta_bl 
  //  << ", original radius: " << radius_unitless
  //  << ", radius in unitless: " << radius_theta_bl_unitless
  //  << ", ratio in unitless: " << radius_theta_bl_unitless/radius_unitless
  //  << std::endl;

  //Step 5: Return this value.
  return radius_theta_bl_unitless * radius_scale;
}

double BaubockOblModel::Dtheta_R( const double& costheta )  {
  // Return dR(theta) / dtheta in "nounits".
  // note that the user supplies cos(theta) and not theta.

  //Step 0: Get m, f, req, and theta.
  double mass_unitless = Mass_nounits();
  double omega_unitless = Omega_nounits();
  double radius_unitless = Req_HT;
  double theta = acos(costheta);

  //Step 1: Turn unitless to cgs.
  double mass_cgs = Units::nounits_to_cgs(mass_unitless, Units::MASS);
  double f_cgs = Units::nounits_to_cgs(omega_unitless, Units::INVTIME)/(2.0*Units::PI);
  double radius_cgs = Units::nounits_to_cgs(radius_unitless, Units::LENGTH);

  //Step 2: Turn cgs to mks.
  double mass_mks = mass_cgs/1000;
  double f_mks = f_cgs;
  double radius_mks = radius_cgs/100;

  //Step 3: Use R_theta_BL(double M_0, double R_0, double f, double theta) to get the radii.
  double drdt = dR_dtheta_BL(mass_mks, radius_mks, f_mks, theta);

  //Step 4: Turn the value to cgs, and then to unitless.
  double drdt_cgs = drdt * 100;
  double drdt_unitless = Units::cgs_to_nounits(drdt_cgs, Units::LENGTH);

  //Step 5: Return this value.
  return drdt_unitless * radius_scale;
}

extern "C" OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model) {
    return new BaubockOblModel(req_nounits, mass_nounits, omega_nounits);
}

extern "C" OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name) {
    return new BaubockOblModel(req_nounits, mass_nounits, omega_nounits);
}
