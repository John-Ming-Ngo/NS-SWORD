/*
Utility class storing assorted widely utilized mathematical equations from the code.
*/

#include "../../include/MathFunctions.h"
#include <cmath>

double P0(const double& mu) {
  return double(1.0);
}

double P2(const double& mu) {
  return double( (3.0*mu*mu - 1.0)/2.0);
}

double P4(const double& mu) {
  return double( (35.0*pow(mu,4.0) - 30.0*mu*mu + 3.0)/8.0 );
}

double Dmu_P0(const double& mu) {
  return double( 0.0 );
}

double Dmu_P2(const double& mu) {
  return double( 3.0*mu ); 
}

double Dmu_P4(const double& mu) {
  return double( mu*(35.0*mu*mu - 15.0)/2.0 );
}
