// Chi.cpp
//
// Computation of Chi-squared

#include "../../include/matpack.h"

#include <exception>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include "../../include/Chi.h"
#include "../../include/OblDeflectionTOA.h"
#include "../../include/OblModelBase.h"
//#include "../../include/PolyOblModelNHQS.h" Unneeded atm.
#include "../../include/PolyOblModelBase.h"
#include "../../include/Exception.h"
#include "../../include/Units.h"
#include "../../include/Struct.h"
#include <stdio.h>
//#define C 2.99792458e5 
using namespace std;


// Given M/R, i and theta, figure out all other angles
class LightCurve ComputeAngles(class LightCurve* incurve, class OblDeflectionTOA* defltoa)			      
{
  class LightCurve curve;
  curve = *incurve;

  double theta, incl, mass, radius, omega, cosgamma;
  double shift_t;
  unsigned int numbins;
  bool infile_is_set;

  theta = curve.para.theta;
  incl = curve.para.incl;
  mass = curve.para.mass;
  radius = curve.para.radius;
  omega = curve.para.omega;
  cosgamma = curve.para.cosgamma;

  curve.eclipse = false;
  curve.ingoing = false;
  curve.problem = false;

  shift_t = curve.junk.shift_t;
  shift_t = 0.001;
  //shift_t = 0.000; //0.001; //before //No apparent change.

  infile_is_set = curve.junk.infile_is_set;
  numbins = curve.numbins;

  double b_guess(0.0);
  double mu, speed;
  double alpha, sinalpha, cosalpha;
  double b, toa_val; 
  double dpsi_db_val;
  // some values which don't affect things for our purposes:
  double dS(100.0), D(1.0);

  //double dmu, dphi; // SMM

  double phi_0(0.0); /* Azimuthal location of spot */
  //phi_0 = -Units::PI; //was in Charlee's code, IDK why needed but testing. Preliminary results: No apparent difference.
 
  bool ingoing(false);

  // Define Vectors

  std::vector< double > phi_e(numbins);
  std::vector< double > psi(numbins);

  //std::vector< double > b(numbins, 0.0);
  std::vector< double > cosdelta(numbins, 0.0);
  //std::vector< double > cosxi(numbins, 0.0);

 // These are calculated in the second loop.
  std::vector< double > dcosalpha_dcospsi(numbins, 0.0);
  std::vector< double > dcosbeta_dcospsi(numbins, 0.0);

  // vectors for 4-point interpolation
  std::vector< double > psi_k(4,0.0);
  std::vector< double > b_k(4,0.0);


  mu = cos(theta);
  speed = omega*radius*sin(theta) / sqrt( 1.0 - 2.0*mass/radius );

  //std::cout << "speed: " << speed << " omega: " << omega << std::endl;

  dS = curve.para.dS;

  // Compute emission time & phase bins & light bending
  for(unsigned int i(0); i < numbins; i++) {
    // SMM: Added an offset of phi_0
    // SMM: Time is normalized to the spin period so 0 < t_e < 1 
    curve.t[i] = i/(1.0*numbins) + shift_t;
    phi_e.at(i) = phi_0 + (2.0*Units::PI) * curve.t[i];
   
    if (cos(incl) == 0.0 && cos(theta) == 0.0)
      psi.at(i) = phi_e.at(i);
    else //Spherical law of cosines
      psi.at(i) = acos(cos(incl)*cos(theta) + sin(incl)*sin(theta)*cos(phi_e.at(i)));
  }

  //  std::cout << "psimax = " << curve.defl.psi_max << std::endl;
  //The integral equation but backwards
  // test for visibility for each value.
  for(unsigned int i(0); i < numbins; i++) {
    int sign(0);
    double bval;
    bool result(false);
    double b1, b2, psi1(0.0), psi2(0.0);
    double xb;

    int j(0), k;

    if (psi.at(i) < curve.defl.psi_max){
      if ( psi.at(i) > curve.defl.psi_b[j] )
        while ( psi.at(i) > curve.defl.psi_b[j] )
          j++;      
      else {
        while (psi.at(i) < curve.defl.psi_b[j])
          j--;
          j++;
      }
      b1 = curve.defl.b_psi[j-1];
      b2 = curve.defl.b_psi[j];
      psi1 = curve.defl.psi_b[j-1];
      psi2 = curve.defl.psi_b[j];
      k = j-2;
      if (j==1) k = 0;
      if (j==3*NN) k = 3*NN - 3;

      for(j=0; j<4; j++){
        b_k.at(j) = curve.defl.b_psi[k+j];
        psi_k.at(j) = curve.defl.psi_b[k+j];
      }
  
      // 4-pt interpolation to find the correct value of b given psi.
      xb = psi.at(i);
      b_guess = (xb-psi_k.at(1))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*b_k.at(0)/
        ((psi_k.at(0)-psi_k.at(1))*(psi_k.at(0)-psi_k.at(2))*(psi_k.at(0)-psi_k.at(3)))
        +(xb-psi_k.at(0))*(xb-psi_k.at(2))*(xb-psi_k.at(3))*b_k.at(1)/
        ((psi_k.at(1)-psi_k.at(0))*(psi_k.at(1)-psi_k.at(2))*(psi_k.at(1)-psi_k.at(3)))
        +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(3))*b_k.at(2)/
        ((psi_k.at(2)-psi_k.at(0))*(psi_k.at(2)-psi_k.at(1))*(psi_k.at(2)-psi_k.at(3)))
        +(xb-psi_k.at(0))*(xb-psi_k.at(1))*(xb-psi_k.at(2))*b_k.at(3)/
        ((psi_k.at(3)-psi_k.at(0))*(psi_k.at(3)-psi_k.at(1))*(psi_k.at(3)-psi_k.at(2)));
    } // End of psi < psi_max
    
    //if (i == 0) std::cout << bval << ",";

    result = defltoa->b_from_psi( fabs(psi.at(i)), mu, bval, sign, curve.defl.b_max, curve.defl.psi_max, 
				  b_guess, fabs(psi.at(i)), b2, fabs(psi.at(i))-psi2, &curve.problem );

    if(result == false) {
      //  std::cout << "i = " << i 
      //	<< ", Not visible at phi_e = " << 180.0*phi_e.at(i)/Units::PI
      //	<< " (deflection not available)." << std::endl;
      curve.visible[i] = false;
      curve.t_o[i] = curve.t[i] ;

      curve.dOmega_s[i] = 0.0;
      curve.eclipse = true;
    }
    else { // there is a solution
      b = bval;
      if(sign < 0 ){
        ingoing = true;
        curve.ingoing = true;
        //	std::cout << "ingoing!"<< std::endl;
      }
      else if(sign > 0){
        ingoing = false;
      }
      else {
        throw(Exception("OblFluxApp.cpp: sign not returned as + or - with success."));
      }
      //sinalpha =  b * sqrt( 1.0 - 2.0*mass/radius ) / radius;
      //cosalpha = sqrt( 1.0 - sinalpha*sinalpha);
      alpha    = asin(b * sqrt( 1.0 - 2.0*mass/radius ) / radius); // alpha can become negative, we cannot simply take positive square roots.

      if(sign < 0) { // alpha is greater than Pi/2
	      alpha = Units::PI - alpha;
      }
      sinalpha = sin(alpha);
      cosalpha = cos(alpha);
      
      curve.alpha[i] = alpha;
      curve.sinalpha[i] = sinalpha;
      curve.cosalpha[i] = cosalpha;

      //Spherical law of cosines
      cosdelta.at(i) = (cos(incl) - cos(theta)*cos(psi.at(i))) / (sin(theta)*sin(psi.at(i))) ;
     
      if ( (cos(theta) < 0) && (cos(incl) < 0 )) 
        cosdelta.at(i) *= -1.0;

      if(sin(psi.at(i)) != 0.0) {
        if (cos(theta)*cos(incl) >= 0) { // Condition needed to fix the fact that sin(gamma) should be negative on the hemisphere facing away from the observer; and positive on the hemisphere facing the observer, according to dR/dTheta.
          curve.cosbeta[i] = cosalpha*cosgamma 
	        + sinalpha*sqrt(1.0 - pow(cosgamma,2.0))*cosdelta.at(i);
        }
        else {
          curve.cosbeta[i] = cosalpha*cosgamma 
	        - sinalpha*sqrt(1.0 - pow(cosgamma,2.0))*cosdelta.at(i);
        }
      }
      else {
        curve.cosbeta[i] = cosalpha*cosgamma;
      }
      //Getting these two values recorded for comparison and documentation purposes.
      curve.beta[i] = acos(curve.cosbeta[i]); // Valid as long as beta is within 0 to pi, which it should be for a NS and for a flat patch with some normal vector n. Will not apply if we can see the back of the patch for some reason.
      curve.sinbeta[i] = sin(curve.beta[i]);
  
      if(curve.cosbeta[i] < 0.0) {
        curve.visible[i] = false;
      }
      else {
        curve.visible[i] = true;
      }

      if(curve.visible[i]) { // visible 
        curve.cosxi[i] = - sinalpha*sin(incl)*sin(phi_e.at(i))/sin(fabs(psi.at(i))); //Spherical law of sines
        curve.eta[i] = sqrt( 1.0 -speed*speed ) / (1.0 - speed*curve.cosxi[i]);

        if(ingoing) {
          //std::cout << "Ingoing b = " << b << std::endl;
          dpsi_db_val = defltoa->dpsi_db_ingoing( b, mu, &curve.problem );
          toa_val = defltoa->toa_ingoing( b, mu, &curve.problem );
        }
        else {
          //std::cout << "Outgoing b = " << b << std::endl;
          dpsi_db_val = defltoa->dpsi_db_outgoing( b, radius, &curve.problem );
          toa_val = defltoa->toa_outgoing( b, radius, &curve.problem );
        }

        //std::cout << "Done computing TOA " << std::endl;
        curve.t_o[i] = curve.t[i] + (omega*toa_val)/(2.0*Units::PI);
            
        dcosalpha_dcospsi.at(i) = fabs(
          sinalpha
          /cosalpha 
          * sqrt(1.0 - 2.0*mass/radius)
          / (radius * sin(fabs(psi.at(i))) * dpsi_db_val)
        );

        //Not yet implemented.
        //dcosbeta_dcospsi.at(i) = fabs(
        //  sinalpha
        //  /cosalpha 
        //  * sqrt(1.0 - 2.0*mass/radius)
        //  / (radius * sin(fabs(psi.at(i))) * dpsi_db_val)
        //);

        // if newtonian
        //dcosalpha_dcospsi.at(i) = 1.0;

        //std::cout << "dcosalpha/dcospsi: " << dcosalpha_dcospsi.at(i) << ", cosbeta: " << curve.cosbeta[i] << std::endl;

        curve.dOmega_s[i] = (dS/(D*D)) 
          * (1.0 / (1.0 - 2.0*mass/radius)) 
          * curve.cosbeta[i] 
          * dcosalpha_dcospsi.at(i);

        if (std::isinf(curve.dOmega_s[i])) {
          std::cerr << "Compute Angles: inf detected." << std::endl;
          std::cout << "i = " << i << " dO[i] = " <<  curve.dOmega_s[i] 
              <<" cosalpha = " << cosalpha << " dpsi/db = " << dpsi_db_val << std::endl;
        }
      } // end visible
      else { // not visible
        //	  phi_o.at(i) = phi_o.at(i-1) + (2.0*Units::PI/numbins);
        curve.t_o[i] = curve.t[i] ;
        curve.dOmega_s[i] = 0.0;
        curve.cosbeta[i] = 0.0;
        curve.eta[i] = 1.0;

      } // end not visible
    }
    curve.b[i] = b;
    curve.psi[i] = psi.at(i);
    curve.cosdelta[i] = cosdelta.at(i);
    curve.dcosalpha_dcospsi[i] = dcosalpha_dcospsi.at(i);
  }
  //  std::cout << "Done Compute Angles " << std::endl;

  return curve;

} // End Compute Angles

class LightCurve ComputeCurve( class LightCurve* angles)
			      		      
{

  class LightCurve curve;

  std::ifstream input;
  std::ofstream output;

  bool ignore_time_delays(false);
  bool infile_is_set;
  double Gamma, aniso, mass, radius;
  int numbins;
  double Gamma1, Gamma2, Gamma3;
  double bbrat;

  double dOmega_mu(0.0); // solid angle for fixed value of mu


  curve = (*angles);
  Gamma = curve.para.Gamma;
  aniso = curve.para.aniso;
  mass = curve.para.mass;
  radius = curve.para.radius;
  infile_is_set = curve.junk.infile_is_set;
  numbins = curve.numbins;
  bbrat = curve.para.bbrat;

  //std::cout << "aniso (compute curve) = " << aniso << "bbrat = " << bbrat << std::endl;

  std::vector< double > newflux(numbins, 0.0); /* rebinned flux */
  std::vector< double > totflux(numbins, 0.0); /* integrated flux */
  std::vector< double > softbb(numbins, 0.0); /* blackbody soft flux */
  std::vector< double > softcm(numbins, 0.0); /* compton soft flux */

  double softbbave(0.0), softcmave(0.0), hardave(0.0);
 
  for (int i(0);i<numbins;i++){
    if (curve.dOmega_s[i] != 0.0){
      dOmega_mu += curve.dOmega_s[i]; 

    }
  }
 
  return curve;

} // end ComputeCurve