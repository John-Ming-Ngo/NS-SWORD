// matpack.h
//
// Declarations for a small subset of the matpack library.
// The matpack copyright notice is printed below.
// I took a subset of the matpack functions in order to not
// to have to link against the whole library, and because
// the library doesn't compile easily on all platforms.
//

// from matpack, i used:
// common.h
// mproots.h
// 
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

#ifndef MATPACK_H
#define MATPACK_H

#include <cfloat>

namespace MATPACK {

  // common.h
  template <class T> inline T CopySign (T x, T y)
    { return (y < 0) ? ((x < 0) ? x : -x) : ((x > 0) ? x : -x); }

  template <class T> inline T MpMax (T x, T y) 
    { return (x>y)?x:y; }
  
  // mproots.h
  double FindZero (double t1, double t2, double (*function)(double), 
		   double tol = DBL_EPSILON);
  
  double FindZeroB (double t1, double t2, double (*function)(double), 
		    double tol = DBL_EPSILON);

  double FindZeroC (double t1, double t2, double b_guess, double psi_guess, double (*function)(double), 
		    double tol = DBL_EPSILON);

  // adaptivesimpson.cpp
  double AdaptiveSimpson (double a, double b, 
			  double (*f)(double x), 
			  double relerr);

} // namespace MATPACK

#endif // MATPACK_H
