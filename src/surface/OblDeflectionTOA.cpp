// OblDeflectionTOA.cpp
//
// This class represents the relationship between
// impact parameter, deflection, and TOA, given 
// an Oblateness model and a Mass
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

#include "../../include/matpack.h"

#include <limits>
#include <exception>
#include <cmath>
#include <iostream>
#include "../../include/OblDeflectionTOA.h"
#include "../../include/OblModelBase.h"
#include "../../include/Exception.h"
#include "../../include/Units.h"
#include "../../include/Struct.h" //Jan 21
// Globals are bad, but there's no easy way to get around it
// for this particular case (need to pass a member function
// to a Matpack routine which does not have a signature to accomodate
// passing in the object pointer).

class LightCurve curve; //Jan 21

const OblDeflectionTOA* OblDeflectionTOA_object;
double OblDeflectionTOA_b_value;
double OblDeflectionTOA_costheta_value;
double OblDeflectionTOA_psi_value;
double OblDeflectionTOA_b_max_value;
double OblDeflectionTOA_psi_max_value;
double OblDeflectionTOA_b_guess;
double OblDeflectionTOA_psi_guess;

double OblDeflectionTOA_psi_integrand_wrapper( double r, bool *prob ){
  double integrand( OblDeflectionTOA_object->psi_integrand( OblDeflectionTOA_b_value, r ));
  if( std::isnan(integrand) ) {
    std::cout << "integrand is nan!" << std::endl;
    std::cout << "r (km) = " << r  << std::endl;
    std::cerr << "OblDeflectionTOA_psi_integrand_wrapper(): returned NaN." << std::endl;
    *prob = true;
    integrand = -7888.0;
    return integrand;
  }
  return integrand;
}


double OblDeflectionTOA_dpsi_db_integrand_wrapper( double r, bool *prob ) {
  double integrand( OblDeflectionTOA_object->dpsi_db_integrand( OblDeflectionTOA_b_value, r ) );
  
  if( std::isnan(integrand) ) {
    throw std::invalid_argument( "received nan value" );
    std::cout << "dpsi_db integrand is nan!" << std::endl;
    std::cout << "r (km) = " << Units::nounits_to_cgs(r, Units::LENGTH)/1.0e5 << std::endl;
    std::cerr << "dpsi_db(): returned NaN." << std::endl;
    *prob = true;
    integrand = -7888.0;
    return integrand;
  }
  
  return integrand;
}

double OblDeflectionTOA_toa_integrand_wrapper( double r, bool *prob ) {
  double integrand( OblDeflectionTOA_object->toa_integrand( OblDeflectionTOA_b_value, r ) );
  if( std::isnan(integrand) ) {
    std::cout << "toa_minus_toa integrand is nan!" << std::endl;
    std::cout << "r (km) = " << Units::nounits_to_cgs(r, Units::LENGTH)/1.0e5 << std::endl;
    std::cerr << "toa_minus_toa(): returned NaN." << std::endl;
    *prob = true;
    integrand = -7888.0;
    return integrand;
  }
  return integrand;
}

double OblDeflectionTOA_toa_integrand_minus_b0_wrapper( double r, bool *prob ) {
  double integrand( OblDeflectionTOA_object->toa_integrand_minus_b0( OblDeflectionTOA_b_value, r ) );
  if( std::isnan(integrand) ) {
    std::cout << "toa_minus_b0 integrand is nan!" << std::endl;
    std::cout << "r (km) = " << Units::nounits_to_cgs(r, Units::LENGTH)/1.0e5 << std::endl;
    std::cerr << "toa_minus_b0(): returned NaN." << std::endl;
    *prob = true;
    integrand = -7888.0;
    return integrand;
  }
  return integrand;
}

double OblDeflectionTOA_rcrit_zero_func_wrapper( double rc ) {
  return double(OblDeflectionTOA_object->rcrit_zero_func( rc, OblDeflectionTOA_b_value ));
}

double OblDeflectionTOA_b_from_psi_ingoing_zero_func_wrapper( double b ) {

  return double(OblDeflectionTOA_object->b_from_psi_ingoing_zero_func( b, 
								       OblDeflectionTOA_costheta_value,
								       OblDeflectionTOA_psi_value)
		);
}

double OblDeflectionTOA_b_from_psi_outgoing_zero_func_wrapper( double b ) {
  // std::cerr << "b = " << Units::nounits_to_cgs(b, Units::LENGTH)/1.0e5 <<std::endl;
  return double(OblDeflectionTOA_object->b_from_psi_outgoing_zero_func( b, 
									OblDeflectionTOA_costheta_value,
									OblDeflectionTOA_psi_value,
									OblDeflectionTOA_b_max_value,
									OblDeflectionTOA_psi_max_value,
									OblDeflectionTOA_b_guess,
									OblDeflectionTOA_psi_guess
									)
		);
}

// End of global pollution.

const double OblDeflectionTOA::INTEGRAL_EPS = 1.0e-7;
const double OblDeflectionTOA::FINDZERO_EPS = 1.0e-6;
const double OblDeflectionTOA::RFINAL_MASS_MULTIPLE = 1.0e6;
//const double OblDeflectionTOA::DIVERGENCE_GUARD = 2.0e-2; // set to 0 to turn off
const double OblDeflectionTOA::DIVERGENCE_GUARD = 0.0;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N = 10000;
const long int OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_INGOING_N = 400;
const double OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_POWER = 4.0;

OblDeflectionTOA::OblDeflectionTOA( OblModelBase* modptr, const double& mass_nounits ) 
  : r_final( RFINAL_MASS_MULTIPLE * mass_nounits ) {
  this->modptr = modptr;
  mass = mass_nounits;
}

double OblDeflectionTOA::bmax_outgoing( const double& rspot ) const{
  // return the maximal value of b for outgoing rays (in funny units!)
  double R=rspot;
  const double M = get_mass();
  double bmax_outgoing( R/sqrt(1.0 - 2.0*M/R) );

  // Below the Schwarzschild photon sphere, the locally tangential ray does
  // not escape to infinity.  The escape cone is bounded by the critical
  // impact parameter 3*sqrt(3)*M.  Stay infinitesimally below that separatrix
  // so the numerical integral never samples the double root at r=3M.
  if (R <= 3.0*M) {
    const double bcrit = 3.0*sqrt(3.0)*M;
    bmax_outgoing = std::nextafter(bcrit, 0.0);
  }
  
  //std::cout << "bmax_outgoing: R = " << R << " bmax = " << bmax_outgoing << std::endl;

  return bmax_outgoing;
}

double OblDeflectionTOA::bmin_ingoing(const double& rspot, const double& cos_theta ) const{
  // return the value of b corresponding to the most ingoing ray allowed
  // (i.e. the one that just grazes the surface and moves toward axis with b=0).
  costheta_check(cos_theta);

  // Need to return the value of b corresponding to dR/dtheta
  // of the surface.
  //double r( modptr->R_at_costheta(cos_theta) );
  double r(rspot);
  double drdth(modptr->Dtheta_R(cos_theta));//(modptr->Dtheta_R(abs(cos_theta))); //TODO: Resolve this at some point (modptr->Dtheta_R(cos_theta)); //(modptr->Dtheta_R(abs(cos_theta))); //Testing what happens if we enforce positivity.
  //std::cout << "bmin_ingoing: r = " << r << " drdth = " << drdth << std::endl;
  double b( r
	    / sqrt( (1.0 - 2.0*get_mass()/r) + pow(drdth / r , 2.0) ) 
	    );



  // NaN check:
  if(std::isnan(b)) {
    std::cerr << "OblDeflectionTOA::bmin_ingoing(): returned NaN." << std::endl;
    //  std::cout << "bmin_ingoing: r = " << r << " drdth = " << drdth << " b = " << b << std::endl;
    b = -7888.0;
    return b;
  }

  return b;
}

bool OblDeflectionTOA::ingoing_allowed( const double& cos_theta ){
  // ingoing rays can be a consideration for locations where dR/dtheta \neq 0
  costheta_check(cos_theta);
  // A photon initially directed inward from inside the photon sphere cannot
  // turn around and escape without re-entering the star.  Only the outward
  // escape cone is physical there.
  return bool( modptr->R_at_costheta(cos_theta) > 3.0*get_mass()
               && modptr->Dtheta_R(cos_theta) != 0.0);
}

double OblDeflectionTOA::rcrit( const double& b, const double& cos_theta, bool *prob ) const{
  double candidate;
  double dummy;

  //  std::cout << "rcrit: " ;

  double r( modptr->R_at_costheta(cos_theta) );

  //std::cout << " r = " << r << " b = " << b << std::endl;

  costheta_check(cos_theta);
  
  if(b > bmax_outgoing(r)
     || b < bmin_ingoing(r, cos_theta)) {
    std::cerr << "OblDeflectionTOA::rcrit() b out-of-range." << std::endl;
    *prob=true;
    dummy=-7888.0;
    return dummy;
  }
  
  if( b == bmax_outgoing(r) ) {
    return double(r);
  }

  OblDeflectionTOA_object = this;
  OblDeflectionTOA_b_value = b;

  // The ingoing turning point is bounded by the Schwarzschild photon
  // sphere and the *local* emission radius.  The former implementation
  // bracketed it between the star's polar and equatorial radii.  That fails
  // for a perfectly valid ingoing ray whose turning point lies below the
  // polar surface, and became visible when the higher-SDIV tabulated surface
  // shifted the limiting impact parameter by a small amount.
  const double lower = 3.0 * get_mass();
  const double upper = r;
  candidate = MATPACK::FindZero(lower,
				upper,
				OblDeflectionTOA_rcrit_zero_func_wrapper);

  if( !std::isfinite(candidate)
      || candidate < lower
      || candidate > upper ) { // this indicates no physical solution
      std::cout << "rcrit = " << candidate << std::endl;
      std::cout << "lower (3M) = " << lower << std::endl;
      std::cout << "local surface radius = " << upper << std::endl;
      std::cerr << "OblDeflectionTOA::rcrit() returned no solution?" << std::endl;
      *prob=true;
      candidate=-7888.0;
      return candidate;
  }


  return candidate;
}

double OblDeflectionTOA::psi_outgoing( const double& b, const double& rspot,
				       const double& b_max, const double& psi_max, bool *prob) const{

  double dummy;

  if(b > b_max
     || b < 0.0 ) {
    std::cerr << "OblDeflectionTOA::psi_outgoing() b out-of-range." << std::endl;
    *prob=true;
    dummy = -7888.0;
    return dummy;
  }

  if (fabs(b-b_max) < 1e-6 )
    {
      return psi_max;
    }
  else{
    // set globals
  
    OblDeflectionTOA_object = this;
    OblDeflectionTOA_b_value = b;

    double psi(0.0);

    if (b!=0.0)  
    {
      double rsurf = rspot;
  
      //      std::cerr << "calculating psi in psi_outgoing, r = " << rsurf << " b = "<< b << std::endl;
      psi = Integration( rsurf,
			      get_rfinal(),
			      OblDeflectionTOA_psi_integrand_wrapper);
    }				
    return psi;
  }
}

// Compute the largest outgoing deflection angle
double OblDeflectionTOA::psi_max_outgoing( const double& b, const double& rspot, bool *prob ) {
  double dummy;
  if(b > bmax_outgoing(rspot)
     || b < 0.0 ) {
    std::cerr << "OblDeflectionTOA::psi_outgoing() b out-of-range." << std::endl;
    *prob=true;
    dummy = -7888.0;
    return dummy;
  }

  // set globals
  OblDeflectionTOA_object = this;
  OblDeflectionTOA_b_value = b;

  double rsurf = rspot;

  //std::cerr << "calculating psi in psi_max_outgoing" << std::endl;

  double psi = Integration( rsurf,
			    get_rfinal(),
			    OblDeflectionTOA_psi_integrand_wrapper);
					
  return psi;
}






double OblDeflectionTOA::psi_ingoing( const double& b, const double& cos_theta, bool *prob ) const{
  costheta_check(cos_theta);
  

  double rspot( modptr->R_at_costheta(cos_theta) );
  double bmax( bmax_outgoing(rspot));
  double dummy;
  //std::cout << "Entering psi_ingoing: cos_theta = " << cos_theta << " b = " << b << " r = " << rspot << std::endl;

  if( (b-bmax) > 1e-15
     || b < bmin_ingoing(rspot, cos_theta)) {
    std::cerr << "OblDeflectionTOA::psi_ingoing() b out-of-range." << std::endl;
    std::cerr << "b = " << b <<std::endl;
    std::cerr << "bmax_outgoing = " << bmax << std::endl;
    std::cerr << "b-bmax = " << b-bmax << std::endl;
    std::cerr << "bmin_ingoing = " << bmin_ingoing(rspot,cos_theta) << std::endl;
    *prob=true;
    dummy = -7888.0;
    return dummy;
  }
  
  // set globals
  OblDeflectionTOA_object = this;
  OblDeflectionTOA_b_value = b;

  double psi_in;

  if ( fabs(b - bmax) <= 1e-10 ){
    psi_in = this->psi_outgoing( bmax, rspot, 100.0, 100.0, &curve.problem );
    //std::cout << "Final Result psi_ingoing: psi_in = " << psi_in << std::endl; 
  }
  else{
  // See psi_outgoing. Use an approximate formula for the integral near
  // rcrit, and the real formula elsewhere
  
    double rcrit = this->rcrit(b, cos_theta, &curve.problem);

    psi_in =  2.0*sqrt( 2.0*(rspot - rcrit)/(rcrit-3.0*get_mass()));

    psi_in += this->psi_outgoing( b, rspot, 100.0, 100.0, &curve.problem );

    //  std::cout << "Final Result psi_ingoing: rcrit = " << rcrit << " psi_in = " << psi_in << std::endl; 
  
  }

  return double( psi_in );
}

bool OblDeflectionTOA::b_from_psi( const double& psi, const double& cos_theta, double& b, int& rdot, 
				   const double& bmax_out, const double& psi_out_max, 
				   const double& b_guess, const double& psi_guess,
				   const double& b2, const double& psi2, bool *prob) {

  // Major Change to b_from_psi
  // a guess for b is given to this routine
  // if the photon is purely outgoing nothing is changed
  // if the photon is initially ingoing a calculation is done

  costheta_check(cos_theta);

  // return bool indicating whether a solution was found.
  // store result in b
  // store -1 in rdot if its an ingoing solution, +1 otherwise.

  // std::cerr << "bfrompsi entered." << std::endl;

  //  double bmax_out( bmax_outgoing(cos_theta) ),
  //psi_out_max( psi_outgoing(bmax_out, cos_theta) ),
  //bcand;

  double bcand;
  double dummyb;
  double dummypsi;

    OblDeflectionTOA_b_max_value = bmax_out;
    OblDeflectionTOA_psi_max_value = psi_out_max;
    OblDeflectionTOA_b_guess = b_guess;
    OblDeflectionTOA_psi_guess = psi_guess;

  if(std::isinf(psi_out_max)) {
    std::cerr << "OblDeflectionTOA::b_from_psi() bfrompsi: inf detected." << std::endl;
    *prob=true;
    dummyb = -7888.0;
    return dummyb;
  }

  // std::cerr << "finished initialiser." << std::endl;

  if(psi < 0) {
    std::cerr << "OblDeflectionTOA::b_from_psi wants psi >= 0!" << std::endl;
    *prob=true;
    dummypsi = -7888.0;
    return dummypsi;
  }
  else if(psi == 0.0) {
    b=0.0;
    rdot = 1;
    return true;
  }
  else if(psi == psi_out_max) {
    b = bmax_out;
    rdot = 1;
    return true;
  }
  else if(psi < psi_out_max) { /* normal outgoing case */
    OblDeflectionTOA_object = this;
    OblDeflectionTOA_costheta_value = cos_theta;
    OblDeflectionTOA_psi_value = psi;

    bcand = b_guess;

    if (bcand > bmax_out){
     
      // std::cout << "bcand - bmax = " << bcand - bmax_out << " setting bcand = bmax - epsilon" << std::endl;

      bcand = bmax_out - 1e-6;  
      
    }
    /*
    bcand = MATPACK::FindZeroC(b_guess, 
			       psi_guess, b2, psi2,
			       OblDeflectionTOA_b_from_psi_outgoing_zero_func_wrapper,
			       OblDeflectionTOA::FINDZERO_EPS);
    */    	
    /*
   bcand = MATPACK::FindZero(b_guess, 
			       b2, 
			       OblDeflectionTOA_b_from_psi_outgoing_zero_func_wrapper,
			       			       OblDeflectionTOA::FINDZERO_EPS);		       
    */

    if( fabs(bcand) <= std::numeric_limits< double >::epsilon()
       ) { // this indicates no soln
      std::cerr << "OblDeflectionTOA::b_from_psi() outgoing returned no solution?" << std::endl;
      *prob=true;
      b = -7888.0;
      return b;
    }
    else {
      b = bcand;
      rdot = 1;
      return true;
    }
  }
  else{ // psi > psi_out_max 
    // Test to see if ingoing photons are allowed

    bool ingoing_allowed(this->ingoing_allowed(cos_theta));
    double bmin_in, psi_in_max;

    //  std::cout << "b_from_psi: psi > psi_out_max" << std::endl;

    if(ingoing_allowed) {

      double rspot = modptr->R_at_costheta(cos_theta);
      //std::cout << "ingoing allowed! rspot = "<< rspot << std::endl;

      bmin_in = bmin_ingoing(rspot, cos_theta);
      //std::cout << "bmin = "<< bmin_in<<std::endl;

      psi_in_max = psi_ingoing(bmin_in, cos_theta, &curve.problem);

      //std::cout << "psimin = "<< psi_in_max<<std::endl;


      if(psi > psi_in_max) {
	return false;
      }
      else if(psi == psi_in_max) {
	b = bmin_in;
	rdot = -1;
	return true;
      }
      else { // psi_out_max < psi < psi_in_max
	OblDeflectionTOA_object = this;
	OblDeflectionTOA_costheta_value = cos_theta;
	OblDeflectionTOA_psi_value = psi;
	//	std::cout << "Compute b for ingoing photon. bmin =  " << bmin_in << " bmax = " << bmax_out << std::endl;
	bcand = MATPACK::FindZero(bmin_in, 
				  bmax_out,
			      OblDeflectionTOA_b_from_psi_ingoing_zero_func_wrapper,
			      OblDeflectionTOA::FINDZERO_EPS);
	//std::cout << "Finished computing b for ingoing photon " << std::endl;
	if( fabs(bmin_in - bcand) <= std::numeric_limits< double >::epsilon()
	    || fabs(bmax_out - bcand) <= std::numeric_limits< double >::epsilon() ) { // this indicates no soln
	  std::cerr << "OblDeflectionTOA::b_from_psi() ingoing returned no solution?" << std::endl;
	  *prob=true;
	  b = -7888.0;
	  return b;
	}
	else {
	  b = bcand;
	  rdot = -1;
	  return true;
	}
      }

    }
    else {
      // std::cerr << "Something is really wrong in b-from-psi!!" << std::endl;
      return false;
    }
  }
  std::cerr << "OblDeflectionTOA::b_from_psi() reached end of function?" << std::endl;
  *prob=true;
  return false;
}


double OblDeflectionTOA::dpsi_db_outgoing( const double& b, const double& rspot, bool *prob ){
  // warning: this integral will diverge near bmax_out. You can't remove
  // the divergence in the same manner as the one for the deflection.
  // check the output!
  //
  //    costheta_check(cos_theta);
  
  double rsurf;
  double dummy;

  if(b > bmax_outgoing(rspot)
     || b < 0.0 ) {
    std::cerr << "OblDeflectionTOA::dpsi_db_outgoing() b out-of-range." << std::endl;
    std::cerr << "b=" << b << "bmax=" << bmax_outgoing(rspot)<< std::endl;
    *prob=true;
    dummy = -7888.0;
    return dummy;
  }
  
  // set globals
  OblDeflectionTOA_object = this;
  OblDeflectionTOA_b_value = b;

  //double rsurf = modptr->R_at_costheta(cos_theta);

  rsurf = rspot;

  double dpsidb = Integration( rsurf,
			       get_rfinal(),
			       OblDeflectionTOA_dpsi_db_integrand_wrapper);
  return dpsidb;
}

double OblDeflectionTOA::dpsi_db_ingoing( const double& b, const double& cos_theta, bool *prob ) {
  // warning: i would have preferred to do this directly,
  // but the integrand diverges at rcrit, and so it is difficult
  // to see how to get the boundary term of this derivative to
  // come out properly. however, we have a closed-form analytical
  // approximation of the integral between rcrit and Rsurf, 
  // so we use it instead for the derivative.
  // 

  //double dumb_val;
  double dummyb;
  costheta_check(cos_theta);
  double rsurf = modptr->R_at_costheta(cos_theta);

  if(b > bmax_outgoing(rsurf)
     || b < bmin_ingoing(rsurf, cos_theta)) {
    std::cerr << "OblDeflectionTOA::dpsi_db_ingoing() b out-of-range." << std::endl;
    std::cerr << "b = " << b <<std::endl;
    std::cerr << "bmax_outgoing = " << bmax_outgoing(rsurf) << std::endl;
    std::cerr << "bmin_ingoing = " << bmin_ingoing(rsurf,cos_theta) << std::endl;

    *prob=true;
    dummyb=-7888.0;
    return dummyb;
  }

  // Don't need this since we're not going to integrate.
  // // set globals
  // OblDeflectionTOA_object = this;
  // OblDeflectionTOA_b_value = b;

  double rcrit = this->rcrit(b, cos_theta, &curve.problem);
 
  double m = this->get_mass();
  double drcrit_db = sqrt(1.0 - 2.0*m/rcrit) 
    / (1.0 - ((b/rcrit)*(m/rcrit)
	      /sqrt(1.0 - 2.0*m/rcrit) 
	      )
       );

  
  double dpsidb_in = -sqrt(2.0)*(drcrit_db)*(rsurf - 3.0*m)
    / (sqrt((rsurf - rcrit)/(rcrit-3.0*m))*pow(rcrit - 3.0*m,2.0)
       );

  return double(dpsidb_in + this->dpsi_db_outgoing( b, rsurf, &curve.problem ));
}

double OblDeflectionTOA::toa_outgoing( const double& b, const double& rspot, bool *prob ){
  // costheta_check(cos_theta);
  double dummy;
  if(b > bmax_outgoing(rspot)
     || b < 0.0 ) {
    std::cerr << "OblDeflectionTOA::toa_outgoing() b out-of-range." << std::endl;
    *prob=true;
    dummy=-7888.0;
    return dummy;
  }
  
  // set globals
  OblDeflectionTOA_object = this;
  OblDeflectionTOA_b_value = b;

  // Note: Use an approximation to the integral near the surface
  // to avoid divergence problems, and then use the real
  // integral afterwards. See astro-ph/07030123 for the formula.

  // The idea is to compute the time for a ray emitted from the
  // surface, and subtract the time for a b=0 ray from r=rpole.
  // To avoid large numbers, part of the subtraction is accomplished
  // in the integrand.

  // Note that for a b=0 ray, Delta T = int(1/(1-2*M/r), r),
  // which is r + 2 M ln (r - 2M)

  //  double rsurf = modptr->R_at_costheta(cos_theta);
 
  double rsurf = rspot;

  double rpole = modptr->R_at_costheta(1.0);

  //std::cout << "toa_outgoing: rpole = " << rpole << std::endl;

  double toa = Integration( rsurf,
			    get_rfinal(),
			    OblDeflectionTOA_toa_integrand_minus_b0_wrapper);

  double toa_b0_polesurf = 
    (rsurf - rpole)
    + 2.0*get_mass()*( log(rsurf - 2.0*get_mass())
		       - log(rpole - 2.0*get_mass())
		       );

  return double( toa - toa_b0_polesurf);
}

double OblDeflectionTOA::toa_ingoing( const double& b, const double& cos_theta, bool *prob ){
  double dummy;
  costheta_check(cos_theta);

  double rsurf = modptr->R_at_costheta(cos_theta);

  if(b > bmax_outgoing(rsurf)
     || b < bmin_ingoing(rsurf,cos_theta)) {
    std::cerr << "OblDeflectionTOA::toa_ingoing() b out-of-range." << std::endl;
    *prob=true;
    dummy=-7888.0;
    return dummy;
  }

  // std::cerr << "DEBUG: rcrit (km) = " << Units::nounits_to_cgs(this->rcrit(b,cos_theta), Units::LENGTH)/1.0e5 << std::endl;

  OblDeflectionTOA_object = this;
  OblDeflectionTOA_b_value = b;

  // See psi_ingoing. Use an approximate formula for the integral near
  // rcrit, and the real formula elsewhere
  
  double rcrit = this->rcrit(b, cos_theta, &curve.problem);
 
  // double toa_in = 2.0 * Integration( rcrit,
  //modptr->R_at_costheta(cos_theta),
  //				     OblDeflectionTOA_toa_integrand_wrapper,
  //				     OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_INGOING_N);

  // double psi_in =  2.0*sqrt( 2.0*(modptr->R_at_costheta(cos_theta) - rcrit)/(rcrit-3.0*get_mass()));

  double toa_in = rcrit/sqrt(1-2.0*get_mass()/rcrit) * 
    2.0*sqrt( 2.0*(modptr->R_at_costheta(cos_theta) - rcrit)/(rcrit-3.0*get_mass()));


//  std::cout << " b = " << b
//	    << " rcrit = " << rcrit
//	    << " toa_in = " << toa_in
//	    << " toa_in = " << rcrit/sqrt(1-2.0*get_mass()/rcrit) * 
//2.0*sqrt( 2.0*(modptr->R_at_costheta(cos_theta) - rcrit)/(rcrit-3.0*get_mass()))
//	    << " toa_out " << this->toa_outgoing(b,cos_theta)
//	    << std::endl;

  double rspot( modptr->R_at_costheta(cos_theta) );

  return double(toa_in + this->toa_outgoing( b, rspot, &curve.problem ));

}


double OblDeflectionTOA::psi_integrand(const double& b, const double& r) const{
  double integrand( b/(r*r*sqrt( 1.0 - pow(b/r,2.0)*(1.0-2.0*get_mass()/r ) ) ) );
  return integrand;
}

double OblDeflectionTOA::dpsi_db_integrand(const double& b, const double& r) const{
  double sqrt_inside = 1.0 - pow(b/r,2.0)*(1.0 - 2.0*get_mass()/r);
  if (sqrt_inside < 0) {
    throw std::invalid_argument( "Square root of negative number occurred; consider b, r:" );
  }
  double integrand( (1.0/(r*r*sqrt(1.0 - pow(b/r,2.0)*(1.0 - 2.0*get_mass()/r))))
		    +
		    ( b*b*(1.0 - 2.0*get_mass()/r) 
		      / (pow(r,4.0)*pow(1.0 - pow(b/r,2.0)*(1.0 - 2.0*get_mass()/r),1.5) ) 
		      )
		    );
  return integrand;
}

double OblDeflectionTOA::toa_integrand(const double& b, const double& r) const{
  double integrand( (1.0/(1.0 - 2.0*get_mass()/r))*((1.0/sqrt( 1.0 - pow(b/r,2.0)*(1.0-2.0*get_mass()/r ) ) ) ) );
  return integrand;
}

double OblDeflectionTOA::toa_integrand_minus_b0(const double& b, const double& r) const{
  double integrand( (1.0/(1.0 - 2.0*get_mass()/r))*((1.0/sqrt( 1.0 - pow(b/r,2.0)*(1.0-2.0*get_mass()/r ) ) ) - 1.0 ) );
  return integrand;
}


double OblDeflectionTOA::rcrit_zero_func( const double& rc, const double& b ) const{
  return double(rc - b*sqrt(1.0 - 2.0*get_mass()/rc ) );
}

double OblDeflectionTOA::b_from_psi_ingoing_zero_func( const double& b, const double& cos_theta, const double& psi ) const{
  return double( psi - this->psi_ingoing(b, cos_theta, &curve.problem) );
}

double OblDeflectionTOA::b_from_psi_outgoing_zero_func( const double& b, const double& cos_theta, const double& psi,
							const double& b_max, const double& psi_max, 
							const double& b_guess, const double& psi_guess) const{
  return double( psi - this->psi_outgoing(b, cos_theta, b_max, psi_max, &curve.problem) );
}

double OblDeflectionTOA::TrapezoidalInteg( const double& a, const double& b, double (*func)(double x, bool *prob),
					   const long int& N ) {
  if( a==b ) return double(0.0);

  double integral(0.0);
  // const long int N(OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_N);
  const double power(OblDeflectionTOA::TRAPEZOIDAL_INTEGRAL_POWER);

  for(int i(1); i <= N; i++) {
    double left, right, mid, width, fmid;
    left = TrapezoidalInteg_pt(a, b, power, N, i-1, &curve.problem);
    right = TrapezoidalInteg_pt(a, b, power, N, i, &curve.problem);
    mid = (left + right)/2.0;
    width = right - left;
    fmid = (*func)(mid, &curve.problem); //Double Check
    
    integral += width*fmid;
  }

  return integral;
}

// power spacing of subdivisions for above integrator, i ranges from 0 to N.
double OblDeflectionTOA::TrapezoidalInteg_pt( const double& a, const double& b, 
					      const double& power, const long int& N,
					      const long int& i, bool *prob) {
  double dummy;
  if(N <= 0){
    std::cerr << "OblDeflectionTOA::TrapezoidalInteg_pt: N <= 0." << std::endl;
    *prob=true;
    dummy=-7888.0;
    return dummy;
  }
  if(i < 0 || i > N){
    std::cerr << "OblDeflectionTOA::TrapezoidalInteg_pt: i out of range." << std::endl;
    *prob=true;
    dummy=-7888.0;
    return dummy;
  }
  return double( a + (b-a)*pow(1.0*i/N,power) );
}

inline double OblDeflectionTOA::Integration( const double& a, const double& b, double (*func)(double x, bool *prob),
					     const long int& N // = TRAPEZOIDAL_INTEGRAL_N
					     ) {
  return TrapezoidalInteg( a, b, func, N );
  // return MATPACK::AdaptiveSimpson( a, b, func, OblDeflectionTOA::INTEGRAL_EPS );
}
