// Chi.cpp
//
// Computation of Chi-squared
#define NDIM 5
#define MPTS 6
#include "PolyOblModelBase.h"

double ChiSquare(class DataStruct* obsdata,
		 class LightCurve* curve,  const double incl, const double aniso, const double theta, const double bbrat, double ts);

class LightCurve ComputeAngles( class LightCurve* incurve,
				class OblDeflectionTOA* defltoa);
			       
class LightCurve ComputeCurve( class LightCurve* angles);


double BlackBody( double T, double E );


int FitCurve(class LightCurve* angles, class DataStruct* obsdata,
	     class OblDeflectionTOA* defltoa, char vert[180], char out_dir[80], double rspot, double omega, double mass, int run, char min_file[180], double ftol);
	      
void amoeba(double p[][NDIM+1], double y[], double ftol, int *nfunk, class LightCurve* curve, class DataStruct* obsdata,
	    class OblDeflectionTOA* defltoa, double omega, double mass, double theta, double rspot, double chi_store[], double i_store[], double a_store[], double t_store[], double b_store[], double p_store[], int *k_value);   //GC

double amotry(double p[][NDIM+1], double y[], double psum[],
	      int ihi, double fac, class LightCurve* curve, class DataStruct* obsdata,
	      class OblDeflectionTOA* defltoa,  double omega, double mass, double theta, double rspot);  //GC -change psum

//class OblDeflectionTOA* recalc(class LightCurve* curve,  double omega, double mass, double theta, double rspot); //Unused atm.
