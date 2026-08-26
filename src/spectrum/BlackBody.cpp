// Routines for Computing BlackBody
// Mostly by Abbie Stevens


#include <exception>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <unistd.h>
#include "../../include/BlackBody.h"
#include "../../include/Units.h"
#include "../../include/Struct.h"
#include "time.h"
#include "../../include/interp.h"
#include <stdio.h>
using namespace std;


// Outputs Specific Intensity vs photon energy for a blackbody spectrum
// Input log_10(T) where T is temperature in Kelvin
// NE is number of photon energies
// Photon energies (in keV) are passed to the function in spectrum.energy
void BBSpectrum(double LogT, class Spectrum* spectrum){

  double Temperature, Energy; // both in keV
  double Intensity;           // erg/cm^2 /(s Hz str)
  int NE = spectrum->Npts;

  Temperature = pow(10,LogT) * Units::K_BOLTZ / Units::EV * 1e-3; // Temperature in keV
  //std::cout << "Temperature: " << Temperature << std::endl;

  //std::cout << "Log(T) = " << LogT
	//    << ";   kT = " << Temperature << " keV"
	//    << std::endl;

  for (unsigned int i(0); i< NE; i++){

    Energy = spectrum->energy[i]; // Energy in keV;

    Intensity = BlackBody(Temperature,Energy);

    //    std::cout << "i = " << i << " energy = " << Energy << std::endl;

    spectrum->intensity[i] = Intensity;
    
  }

  
}




/**************************************************************************************/
/* Blackbody:                                                                         */
/*           computes the monochromatic blackbody flux in units of erg/cm^2			  */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E = monochromatic energy in keV * redshift / eta                             */
/**************************************************************************************/
double BlackBody( double T, double E ) {   // Blackbody flux in units of erg/cm^2
  return ( 2.0e9 / pow(Units::C * Units::H_PLANCK, 2) * pow(E * Units::EV, 3) / (exp(E/T) - 1) ); 

  // return (pow(E/T,3.0)/(exp(E/T) - 1));

    // the e9 is to switch E from keV to eV; Units::EV gets it from eV to erg, since it's first computed in erg units.
    // the switch from erg units to photon count units happens above just after this is called.
} // end Blackbody






/**************************************************************************************/
/* EnergyBandFlux:                                                                    */
/*                computes the blackbody flux in units of erg/cm^2 using trapezoidal  */
/*                rule for approximating an integral                                  */
/*                variant of Bradt equation 6.6                                       */
/*                T, E1, E2 put into eV in this routine                               */
/*																					  */
/* pass: T = the temperature of the hot spot, in keV                                  */
/*       E1 = lower bound of energy band in keV * redshift / eta                      */
/*       E2 = upper bound of energy band in keV * redshift / eta                      */
/**************************************************************************************/
double EnergyBandFlux( double T, double E1, double E2 ) {
	T *= 1e3; // from keV to eV
	// x = E / T
	E1 *= 1e3; // from keV to eV
	E2 *= 1e3; // from keV to eV
	
	/********************************************/
   	/* VARIABLE DECLARATIONS FOR EnergyBandFlux */
    /********************************************/
	
	// a, b, x, n, h as defined by Mathematical Handbook eqn 15.16 (Trapezoidal rule to approximate definite integrals)
	double a = E1 / T;          // lower bound of integration
	double b = E2 / T;          // upper bound of integration
	double current_x(0.0);      // current value of x, at which we are evaluating the integrand; x = E / T; unitless
	unsigned int current_n(0);  // current step
	//unsigned int n_steps(900); // total number of steps
	unsigned int n_steps(800);
	// This number of steps (100) is optimized for Delta(E) = 0.3 keV
	double h = (b - a) / n_steps;     // step amount for numerical integration; the size of each step
	double integral_constants = 2.0 * pow(T*Units::EV,3) / pow(Units::C,2) / pow(Units::H_PLANCK,3); // what comes before the integral when calculating flux using Bradt eqn 6.6 (in units of photons/cm^2/s)
	double flux(0.0);           // the resultant energy flux density; Bradt eqn 6.17
	
	// begin trapezoidal rule
	current_x = a + h * current_n;
	flux = Bradt_flux_integrand(current_x);

	for ( current_n = 1; current_n < n_steps-1; current_n++ ) {
		current_x = a + h * current_n;
		flux += 2.0 * Bradt_flux_integrand(current_x);
	}
	
	current_x = a + h * current_n;
	flux += Bradt_flux_integrand(current_x);
	
	flux *= h/2.0;	
	// end trapezoidal rule; numerical integration complete!
	
	flux *= integral_constants;

	return flux;
} // end EnergyBandFlux



/**************************************************************************************/
/* Bradt_flux_integrand:                                                              */
/*                      integrand of Bradt eqn 6.6 when integrating over nu, modified */
/*                      so the exponent is 2 not 3, so that it comes out as photon    */
/*                      number flux, instead of erg flux                              */
/*																					  */
/* pass: x = current_x from above routine                                             */
/**************************************************************************************/
double Bradt_flux_integrand( double x ) {
	return ( pow(x,2) / (exp(x) - 1) );  // 2 (not 3) for photon number flux
} // end Bradt_flux_integrand

double Hopf( double cosalpha) {

  return (0.42822+0.92236*cosalpha-0.085751*pow(cosalpha,2)) ;

}
/*
Stephan-Boltzmann law emissivity
Returns bolometric flux in ergs/cm^2
*/
double StephBoltzmann(double T_k, double emissivity) {
	return Units::STEPH_BOLTZMANN * emissivity * pow(T_k, 4);
}
