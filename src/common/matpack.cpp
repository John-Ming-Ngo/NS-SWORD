// matpack.h
//
// Source code for a small subset of the matpack library.
// The matpack copyright notice is printed below.
// I took a subset of the matpack functions in order to not
// to have to link against the whole library, and because
// the library doesn't compile easily on all platforms.
//

// from matpack, i used:
// common.h
// mproots.h
// adaptsimpson.cpp
// zero.cpp
// zerob.cpp

/*-----------------------------------------------------------------------------*\
| common definitions include file                                      common.h |
|                                                                               |
| Last change: Sep 12, 2004                                                     |
|                                                                               |
| Matpack Library Release 1.8.0                                                 |
| Copyright (C) 1991-2004 by Berndt M. Gammel. All rights reserved.             |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.                                                |
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.                                                    |
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without prior consent of the author.                       |
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration |
| and installation of Matpack.                                                  |
|                                                                               |
\*-----------------------------------------------------------------------------*/

#include "../../include/matpack.h"
#include "../../include/Exception.h"
#include <cmath>
#include <cfloat>
#include <iostream>
/*-----------------------------------------------------------------------------*\
| Matpack univariate zero finder                                       zero.cpp |
|                                                                               |
| Last change: Sep 12, 2004                                                     |
|                                                                               |
| Matpack Library Release 1.8.0                                                 |
| Copyright (C) 1991-2004 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its  entirety and its | 
| documentation for  non-commercial purpose and  without fee is hereby granted, | 
| provided that this license information and copyright notice appear unmodified | 
| in all copies. This software is provided 'as is'  without  express or implied | 
| warranty.  In no event will the author be held liable for any damages arising | 
| from the use of this software.						|
| Note that distributing Matpack 'bundled' in with any product is considered to | 
| be a 'commercial purpose'.							|
| The software may be modified for your own purposes, but modified versions may | 
| not be distributed without prior consent of the author.			|
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration	|
| and installation of Matpack.							|
|                                                                               |
\*-----------------------------------------------------------------------------*/

// #include "mproots.h"

namespace MATPACK {

//----------------------------------------------------------------------------//
//
// double FindZero (double t1, double t2, 
//                  double (*function)(double), 
//                  double tol = DBL_EPSILON);
//
// Univariate zero finding. Computes the zero of a given function f(t) 
// between t1 and t2 to a given tolerance.
//
// Arguments:
//
//   double t1			lower boundary of interval
//   double t2			upper boundary of interval
//   double (*function)(double)	the given function 
//   double tol			acceptable tolerance for root position, if 
//                              omitted then the machine precision DBL_EPSILON 
//                              is used
//
// Method:
//
//   The method is a combination of the regula falsi and the midpoint method.
//   It is a modified version of the VIM - (control data user group) routine
//   with catalog identification "c2bkyzero" written by Loren P. Meissner, 1965.
//
// Notes:
//   
//   1) It doesn't matter if t1 is smaller or larger than t2.
//   2) If there is no zero within t1, t2 then a boundary value 
//      (within the tolerance) is returned.
//   3) This algorithm is more than 1.6 as fast as Brent's zero finder.
//      For a simple problem like the double root of a 4th order polynomial.
//   4) The algorithm has been completely rewritten in C++ by B.M. Gammel, 1995.
//
//----------------------------------------------------------------------------//

double FindZero (double t1, double t2, double (*function)(double), double tol)
{
    double a, b, c, s, fa, fb, fc, fs, e, g, h, y,  fg, fy;

    /*
    if (tol <= 0) Matpack.Error(Mat::ArgumentDomain,
				"FindZero: Tolerance must be positve");
    */
    // Local modification:
    if (tol <= 0) throw(Exception("FindZero: Tolerance must be positive"));


    a = t1;
    b = t2;

    fa = function(a);
    fb = function(b);

   

    if (std::isnan(fa) || std::isnan(fb) ){
      std::cout << "Matpack: FindZero: a = " << a << " f(a) = " << fa
		<< " b = " << b << " f(b) = " << fb << std::endl;
      throw(Exception("FindZero: nan encountered!"));
     
    }

    c = a;
    fc = fa;
    s = c;
    fs = fc;

    for (;;) { // main loop 

	h = (b + c) * 0.5;
	if (fabs(h - b) <= tol) return h;
	
	if (fabs(fb) <= fabs(fc)) {
	    y = s;
	    fy = fs;
	    g = c;
	    fg = fc;
	    s = b;
	    fs = fb;
	} else { 
	    y = b;
	    fy = fb;
	    g = b;
	    fg = fb;
	    s = c;
	    fs = fc;
	}
	
	if (fy != fs) {    
	    e = (s * fy - y * fs) / (fy - fs);
	    if (fabs(e-s) <= tol) e = s + CopySign(tol,g-s);
	    b = ((e - h) * (s - e) < 0.0) ? h : e;
	} else
	    b = h;

	// call function
	fb = function(b);
	
	if (fg * fb < 0.0) {  
	    c = g;
	    fc = fg;
	} else {
	    c = s;
	    fc = fs;
	}
    } // end of main loop
}

} // namespace MATPACK

//----------------------------------------------------------------------------//



/*-----------------------------------------------------------------------------*\
| Matpack univariate zero finder                                      zerob.cpp |
|                                                                               |
| Last change: Sep 12, 2004                                                     |
|                                                                               |
| Matpack Library Release 1.8.0                                                 |
| Copyright (C) 1991-2004 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its  entirety and its | 
| documentation for  non-commercial purpose and  without fee is hereby granted, | 
| provided that this license information and copyright notice appear unmodified | 
| in all copies. This software is provided 'as is'  without  express or implied | 
| warranty.  In no event will the author be held liable for any damages arising | 
| from the use of this software.						|
| Note that distributing Matpack 'bundled' in with any product is considered to | 
| be a 'commercial purpose'.							|
| The software may be modified for your own purposes, but modified versions may | 
| not be distributed without prior consent of the author.			|
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration	|
| and installation of Matpack.							|
|                                                                               |
\*-----------------------------------------------------------------------------*/

// #include "mproots.h"

namespace MATPACK {

//----------------------------------------------------------------------------//
//
// double FindZeroB (double t1, double t2, 
//                   double (*function)(double t), 
//                   double tol = DBL_EPSILON);
//
// Univariate zero finding. Computes the zero of a given function f(t) 
// between t1 and t2 to a given tolerance.
//
// Arguments:
//
//   double t1			  lower boundary of interval
//   double t2			  upper boundary of interval
//   double (*function)(double t) the given function 
//   double tol			  acceptable tolerance for root position, if 
//                                omitted then the machine precision DBL_EPSILON 
//                                is used
// Method:
//
//   This is "Brent's zero finder", a slightly  modified  translation  of
//   the Algol 60 procedure zero given in:
//     Richard Brent, "Algorithms for
//     minimization without derivatives", Prentice - Hall, Inc. (1973).
//   See also:
//     G.Forsythe, M.Malcolm, C.Moler, Computer methods for mathematical
//     computations. M., Mir, 1980, p.180 of the Russian edition
//
//   The function makes use of the bissection procedure combined with
//   the linear or quadratic inverse interpolation.
//
//   This function returns an approximate location for the root with accuracy
//   4*DBL_EPSILON*abs(x) + tol.
//
//   At every step program operates three abscissae - a, b, and c.
//
//	b - the last and the best approximation to the root
//	a - the last but one approximation
//	c - the last but one or even earlier approximation such that
//		1) |f(b)| <= |f(c)|
//		2) f(b) and f(c) have opposite signs, i.e. b and c encompass
//		   the root
//
// At every step this function computes two new approximations, one by the 
// bissection procedure and the other one from interpolation (if a,b, and c
// are all different the quadratic interpolation is used, linear otherwise).
// If the latter (i.e. obtained by the interpolation) point looks
// reasonable (i.e. falls within the current interval [b,c] not close
// to the end points of the interval), the point is accepted as a new
// approximation to the root. Otherwise, the bissection result is used.
// Therefore, the range of uncertainty is guaranteed to be reduced at 
// least by the factor of 1.6.
//
//----------------------------------------------------------------------------//

  double FindZeroC (double b1, double p1, double b2, double p2, 
		    double (*function)(double), double tol)
{
    double a, b, c, s, fa, fb, fc, fs, e, g, h, y,  fg, fy;

  

    /*
    if (tol <= 0) Matpack.Error(Mat::ArgumentDomain,
				"FindZero: Tolerance must be positve");
    */
    // Local modification:
    if (tol <= 0) throw(Exception("FindZero: Tolerance must be positive"));


    a = b1;
    b = b2;

    fa = p1;
    fb = p2;

    /*
    fa = function(a);
    fb = function(b);

    
    if (b_guess != 0.0 || b_guess < b){
      z = b_guess;
      fz = psi_guess;
      //fz = function(z);

      if ( fa*fz < 0.0 ){
	b = z;
	fb = fz;
      }
      else{
	a = z;
	fa = fz;
      }
    }
    
    */        


    c = a;
    fc = fa;
    s = c;
    fs = fc;

    for (;;) { // main loop 

	h = (b + c) * 0.5;
	if (fabs(h - b) <= tol) return h;
	
	if (fabs(fb) <= fabs(fc)) {
	    y = s;
	    fy = fs;
	    g = c;
	    fg = fc;
	    s = b;
	    fs = fb;
	} else { 
	    y = b;
	    fy = fb;
	    g = b;
	    fg = fb;
	    s = c;
	    fs = fc;
	}
	
	if (fy != fs) {    
	    e = (s * fy - y * fs) / (fy - fs);
	    if (fabs(e-s) <= tol) e = s + CopySign(tol,g-s);
	    b = ((e - h) * (s - e) < 0.0) ? h : e;
	} else
	    b = h;

	// call function
	fb = function(b);

	//	std::cout << "b=" << b << " fb= " << fb << std::endl;
	
	if (fg * fb < 0.0) {  
	    c = g;
	    fc = fg;
	} else {
	    c = s;
	    fc = fs;
	}
    } // end of main loop
}

 // namespace MATPACK



double FindZeroB (double t1, double t2, double (*function)(double), double tol)
{
    double a, b, c, fa, fb, fc;	

    /*
    if (tol <= 0) Matpack.Error(Mat::ArgumentDomain,
				"FindZeroB: Tolerance must be positive");
    */
    // Local modification:
    if (tol <= 0) throw(Exception("FindZeroB: Tolerance must be positive"));

    if (t1 < t2) {
	a = t1; 
	b = t2;
    } else if (t1 == t2) {
	a = t1;
	b = t2*(1+DBL_EPSILON);
    } else {
	a = t2;
	b = t1;
    }

    fa = function(a);  
    fb = function(b);

    c = a;   
    fc = fa;

    for (;;) {				// Main iteration loop
  
	double prev_step = b-a;		// Step from the previous iteration
   
	if (fabs(fc) < fabs(fb)) {
	    a = b;  b = c;  c = a;          // Swap data so that b would be the
	    fa=fb;  fb=fc;  fc=fa;	// best approximation found so far
	}
					// Estimate the effective tolerance
	double tol_act = 2*DBL_EPSILON*fabs(b) + 0.5*tol;
	double new_step = (c-b)/2;	// Bissection step for this iteration

	if (fabs(new_step) <= tol_act || fb == 0)
	    return b;			// Acceptable approximation is found

					// Figuring out if the interpolation
					// can be tried
	if (fabs(prev_step) >= tol_act	// If prev_step was large enough
	     && fabs(fa) > fabs(fb)) {	// and was in true direction,
					// Interpolatiom may be tried

	    double p;      		// Interpolation step is calcu-
	    double q;      		// lated in the form p/q; divi-
					// sion operations is delayed
					// until the last moment
	    double cb = c-b;

	    if (a == c) {		// If we've got only two distinct
					// points linear interpolation
		double s1 = fb/fa;	// can only be applied
		p = cb*s1;
		q = 1.0 - s1;

	    } else {			// Quadratic inverse interpolation

		double s1, s2;
		q = fa/fc;  
		s1 = fb/fc;  
		s2 = fb/fa;
		p = s2 * (cb*q*(q-s1) - (b-a)*(s1-1.0));
		q = (q-1.0) * (s1-1.0) * (s2-1.0);
	    }

	    if (p > 0)			// Formulas above computed new_step
		q = -q;			// = p/q with wrong sign (on purpose).
	    else			// Correct this, but in such a way so
		p = -p;			// that p would be positive
      
	    if (p < (0.75*cb*q-fabs(tol_act*q)/2) // If b+p/q falls in [b,c]
		&& p < fabs(prev_step*q/2))       // and isn't too large
		new_step = p/q;		// it is accepted
					// If p/q is too large then the
					// bissection procedure can
					// reduce [b,c] to a larger extent
	}

	if (fabs(new_step) < tol_act)	// Adjust the step to be not less
	    if ( new_step > 0)		// than the tolerance
		new_step = tol_act;
	    else
		new_step = -tol_act;

	a = b;  			// Save the previous approximation
	fa = fb;	
	b += new_step;  		// Do step to a new approximation
	fb = function(b);
	if ((fb > 0 && fc > 0) || (fb < 0 && fc < 0)) {
					// Adjust c for it to have the sign
	    c = a;  fc = fa;		// opposite to that of b
	}
    }
}

} // namespace MATPACK

//----------------------------------------------------------------------------//

/*-----------------------------------------------------------------------------*\
| Matpack integration package - adaptive simpson               adaptsimpson.cpp |
|                                                                               |
| Last change: Sep 12, 2004                                                     |
|                                                                               |
| Matpack Library Release 1.8.0                                                 |
| Copyright (C) 1991-2004 by Berndt M. Gammel                                   |
|                                                                               |
| Permission to  use, copy, and  distribute  Matpack  in  its entirety  and its |
| documentation  for non-commercial purpose and  without fee is hereby granted, |
| provided that this license information and copyright notice appear unmodified |
| in all copies.  This software is provided 'as is'  without express or implied |
| warranty.  In no event will the author be held liable for any damages arising |
| from the use of this software.                                                |
| Note that distributing Matpack 'bundled' in with any product is considered to |
| be a 'commercial purpose'.                                                    |
| The software may be modified for your own purposes, but modified versions may |
| not be distributed without prior consent of the author.                       |
|                                                                               |
| Read the  COPYRIGHT and  README files in this distribution about registration |
| and installation of Matpack.                                                  |
|                                                                               |
\*-----------------------------------------------------------------------------*/

// #include "common.h"

namespace MATPACK {

//-----------------------------------------------------------------------------//
// double AdaptiveSimpson (double a, double b, 
//                         double (*f)(double x), 
//                         double relerr);
//-----------------------------------------------------------------------------//
// 
// DESCRIPTION:
//
//   This routine integrates the function f in the intervall [a,b] to a relative
//   accuracy of relerr using an adaptive Simpson algorithm. This is a general
//   purpose algorithm which can be applied to a variety of problems. It is
//   applicable also for functions that have singularities or which are rapidly
//   changing. But it is not very fast. The result becomes exact for polynomials
//   of degree three.
//
//   The error boundary should be considered carefully, because two 
//   approximations are compared. If there is a strong cancelation of terms 
//   of the integral, for instance if the sign of the function is oscillating 
//   rapidely, then the relative error might be larger than expected.
//
//   In principal the adaptive quadrature should be used only in the case
//   of complicated integrals (for instance with singlarities or fast changing
//   functions) or as a test of the results of faster specialized integration
//   methods (Gauss integration or other non-adaptive methods).
//
//   The function  is based on algorithm 145 of  W. M. McKeeman
//   from the "Collected Algorithm from ACM". 
//
// IMPLEMENTATION NOTES:
//
//   Recursion is explicitly used, since it is quite effective
//   in the C language. A non-recursive Fortran version can be found in
//   ACM algorithm 182.  
//  
//   C++ version by Berndt M. Gammel, last change Nov 13, 1996.
//
//-----------------------------------------------------------------------------//

const int maxlev = 20;  // maximum recursion level - should be sufficient always

//-----------------------------------------------------------------------------//
// local function prototype
//-----------------------------------------------------------------------------//

static double RecursiveSimpson (double a, double da, 
                                double fa, double fm, double fb,
                                double area, double est, double relerr,
                                double (*funct)(double x),
                                int &level, int &levmax);

//-----------------------------------------------------------------------------//

double AdaptiveSimpson (double a, double b, 
                        double (*funct)(double x), 
                        double relerr)
{
    int levmax = 1, level = 1;
    
    return RecursiveSimpson(a, 
                            b - a, 
                            funct(a),
                            funct(0.5*(a + b)) * 4, 
                            funct(b), 
                            1.0, 1.0, fabs(relerr), funct, level, levmax);
}


//-----------------------------------------------------------------------------//

static double RecursiveSimpson (double a, double da, 
                                double fa, double fm, double fb,
                                double area, double est, double relerr,
                                double (*funct)(double x),
                                int &level, int &levmax)
{
    const double norm = 0.588,          // heuristic constant = 1/sqrt(3)
                 one_third = 1.0 / 3.0,
                 one_18th  = 1.0 / 18.0;

    double dx  = one_third * da,
           ddx = one_18th * da,
           relerrs = relerr * norm,
           x1 = a + dx,
           x2 = x1 + dx,
           f1 = 4 * funct(a + 0.5 * dx),
           f2 = funct(x1),
           f3 = funct(x2),
           f4 = 4 * funct(a + 2.5 * dx),
           est1 = (fa + f1 + f2) * ddx,
           est2 = (f2 + fm + f3) * ddx,
           est3 = (f3 + f4 + fb) * ddx,
           absarea = area - fabs(est) + fabs(est1) + fabs(est2) + fabs(est3),
           sum = est1 + est2 + est3;

    level++;
    levmax = MpMax(level,levmax);
    if ( ((fabs(est - sum) > relerr*absarea) || (est == 1.0)) && (level < maxlev) )
      sum = RecursiveSimpson( a,dx,fa,f1,f2,absarea,est1,relerrs,funct,level,levmax)
          + RecursiveSimpson(x1,dx,f2,fm,f3,absarea,est2,relerrs,funct,level,levmax)
          + RecursiveSimpson(x2,dx,f3,f4,fb,absarea,est3,relerrs,funct,level,levmax);
    level--;
    return sum;
}

} // namespace MATPACK

//-----------------------------------------------------------------------------//
