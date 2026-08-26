// Struct.h


#ifndef STRUCT_H
#define STRUCT_H

#include <exception>
#include <vector>
#include <float.h>

#define NN 20
#define NUMBINS 256
#define NCURVES 2

struct Parameters {
  double theta;
  double dS;
  double incl;
  double aniso;
  double Gamma;
  double mass;
  double radius;
  double omega;
  double cosgamma;
  double bbrat;
  double ts;
};

struct Junk {
  double shift_t;
  bool infile_is_set;
};


class Defl {
 public:
  double psi_b[3*NN+1];
  double b_psi[3*NN+1];
  double psi_max;
  double b_max;
  //class OblDeflectionTOA defltoa;
};

class LightCurve {
 public:
  double t[NUMBINS+1];
  double f[NCURVES+1][NUMBINS+1];
  bool visible[NUMBINS+1];
  double t_o[NUMBINS+1];
  double alpha[NUMBINS+1];
  double sinalpha[NUMBINS+1];
  double cosalpha[NUMBINS+1];

  double beta[NUMBINS+1];
  double sinbeta[NUMBINS+1];
  double cosbeta[NUMBINS+1];
  
  double eta[NUMBINS+1];
  double cosxi[NUMBINS+1];
  double dOmega_s[NUMBINS+1];
  double b[NUMBINS+1];
  double psi[NUMBINS+1];
  double cosdelta[NUMBINS+1];
  double dcosalpha_dcospsi[NUMBINS+1];
  struct Parameters para;
  struct Junk junk;
  class Defl defl;
  unsigned int numbins;
  bool eclipse; // True if an eclipse occurs
  bool ingoing; // True if one or more points are ingoing
  bool problem; // True if a problem occurs
  bool para_read_in; //True if reading in initial parameter values from a file
};


class DataStruct {
 public:
  double t[NUMBINS+1];
  double f[NCURVES+1][NUMBINS+1];
  double err[NCURVES+1][NUMBINS+1];
  double chisquare[NUMBINS+1];
  double shift;
  unsigned int numbins;
};

struct Flags {
  unsigned int spectral_model;  // stating which model we're using -- definitions of models given elsewhere
  bool kelvin; // true means temperature in Kelvin; false means in keV
  bool logEflag; // true means energy bands evenly spaced in log space; false means linear space
};

class Atmo {                     // Stores all the data about the atmosphere!
	public:
	double *nsxinte;  			 // intensity values
	double *nsxangl; // cos(alpha) where alpha = emission angle from normal
	double *nsxloget;                        // log(photonenergy/temperature)
	double *nsxlogTeff;
	double *nsxlogg; // log(surface gravity)
	int Nmu;
	int NlogTeff;
	int Nlogg;
	int NlogE;
	int Npts;

};

class Spectrum {
 public:
  int Npts; // Number of energy values = Atmo.NlogE
  double *energy; // photon enegy values
  double *intensity; // specific intensity for given photon energy
};


#endif // STRUCT_H
