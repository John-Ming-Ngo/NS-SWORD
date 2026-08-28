/***************************************************************************************/
/*                                     Hydrogen.cpp

    This holds the routines used to load and access the Hydrogen Atmosphere Model
    computed by Wynn Ho. 

*/
/********************************************************************************NS*******/
#include <exception>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <sys/stat.h>
#include <unistd.h>
#ifdef _WIN32
#include <direct.h>
#endif
#include "../../include/Hydrogen.h"
#include "../../include/Units.h"
#include "../../include/Struct.h"
#include "time.h"
#include "../../include/interp.h"
#include "../../include/nrutil.h"
#include <stdio.h>
using namespace std;

namespace {

const char* const NSX_H_DEFAULT_SOURCE_PATH = "inputs/atmosphere/nsx_H_v200804.out";
const char* const NSX_H_CACHE_PATH = "build_objs/cache/nsx_H_v200804.cache";
const char NSX_H_CACHE_MAGIC[8] = {'N', 'S', 'X', 'H', 'C', '0', '2', '\0'};

struct NSXSourceIdentity {
  std::uint64_t size;
  std::int64_t mtime;
};

const char* nsx_h_source_path() {
  const char* configured = std::getenv("NS_SWORD_ATMOSPHERE");
  return configured && configured[0] ? configured : NSX_H_DEFAULT_SOURCE_PATH;
}

bool source_identity(const char* path, NSXSourceIdentity* identity) {
  struct stat info;
  if (stat(path, &info) != 0) return false;
  identity->size = static_cast<std::uint64_t>(info.st_size);
  identity->mtime = static_cast<std::int64_t>(info.st_mtime);
  return true;
}

template <typename T>
bool read_binary(std::ifstream& in, T* value) {
  in.read(reinterpret_cast<char*>(value), sizeof(T));
  return static_cast<bool>(in);
}

template <typename T>
bool write_binary(std::ofstream& out, const T& value) {
  out.write(reinterpret_cast<const char*>(&value), sizeof(T));
  return static_cast<bool>(out);
}

std::uint64_t fnv1a_bytes(std::uint64_t hash, const void* data, std::size_t size) {
  const unsigned char* bytes = static_cast<const unsigned char*>(data);
  const std::uint64_t prime = UINT64_C(1099511628211);
  for (std::size_t i = 0; i < size; ++i) {
    hash ^= bytes[i];
    hash *= prime;
  }
  return hash;
}

std::uint64_t atmosphere_checksum(const Atmo* atmo) {
  std::uint64_t hash = UINT64_C(14695981039346656037);
  hash = fnv1a_bytes(hash, atmo->nsxangl,
                     static_cast<std::size_t>(atmo->Nmu) * sizeof(double));
  hash = fnv1a_bytes(hash, atmo->nsxinte,
                     static_cast<std::size_t>(atmo->Npts) * atmo->Nmu * sizeof(double));
  return hash;
}

bool load_nsx_cache(Atmo* atmo, const NSXSourceIdentity& source) {
  std::ifstream in(NSX_H_CACHE_PATH, std::ios::in | std::ios::binary);
  if (!in) return false;

  char magic[8] = {};
  std::uint32_t n_log_teff = 0, n_log_g = 0, n_log_e = 0, n_mu = 0;
  std::uint64_t source_size = 0;
  std::int64_t source_mtime = 0;
  std::uint64_t expected_checksum = 0;
  in.read(magic, sizeof(magic));
  if (!in || std::memcmp(magic, NSX_H_CACHE_MAGIC, sizeof(magic)) != 0 ||
      !read_binary(in, &n_log_teff) || !read_binary(in, &n_log_g) ||
      !read_binary(in, &n_log_e) || !read_binary(in, &n_mu) ||
      !read_binary(in, &source_size) || !read_binary(in, &source_mtime) ||
      !read_binary(in, &expected_checksum)) {
    return false;
  }

  if (n_log_teff != 35 || n_log_g != 14 || n_log_e != 166 || n_mu != 67 ||
      source_size != source.size || source_mtime != source.mtime) {
    return false;
  }

  const std::uint64_t intensity_count =
      static_cast<std::uint64_t>(n_log_teff) * n_log_g * n_log_e * n_mu;
  in.read(reinterpret_cast<char*>(atmo->nsxangl),
          static_cast<std::streamsize>(n_mu * sizeof(double)));
  in.read(reinterpret_cast<char*>(atmo->nsxinte),
          static_cast<std::streamsize>(intensity_count * sizeof(double)));
  return static_cast<bool>(in) && atmosphere_checksum(atmo) == expected_checksum;
}

void ensure_cache_directory() {
#ifdef _WIN32
  if (_mkdir("build_objs/cache") != 0 && errno != EEXIST) return;
#else
  if (mkdir("build_objs/cache", 0755) != 0 && errno != EEXIST) return;
#endif
}

void write_nsx_cache(const Atmo* atmo, const NSXSourceIdentity& source) {
  ensure_cache_directory();
  std::ostringstream temporary_name;
  temporary_name << NSX_H_CACHE_PATH << ".tmp." << static_cast<long>(getpid());
  const std::string temporary_path = temporary_name.str();
  std::ofstream out(temporary_path.c_str(), std::ios::out | std::ios::binary | std::ios::trunc);
  if (!out) return;

  const std::uint32_t n_log_teff = static_cast<std::uint32_t>(atmo->NlogTeff);
  const std::uint32_t n_log_g = static_cast<std::uint32_t>(atmo->Nlogg);
  const std::uint32_t n_log_e = static_cast<std::uint32_t>(atmo->NlogE);
  const std::uint32_t n_mu = static_cast<std::uint32_t>(atmo->Nmu);
  const std::uint64_t intensity_count =
      static_cast<std::uint64_t>(atmo->Npts) * atmo->Nmu;
  const std::uint64_t checksum = atmosphere_checksum(atmo);

  out.write(NSX_H_CACHE_MAGIC, sizeof(NSX_H_CACHE_MAGIC));
  if (!write_binary(out, n_log_teff) || !write_binary(out, n_log_g) ||
      !write_binary(out, n_log_e) || !write_binary(out, n_mu) ||
      !write_binary(out, source.size) || !write_binary(out, source.mtime) ||
      !write_binary(out, checksum)) {
    out.close();
    std::remove(temporary_path.c_str());
    return;
  }
  out.write(reinterpret_cast<const char*>(atmo->nsxangl),
            static_cast<std::streamsize>(atmo->Nmu * sizeof(double)));
  out.write(reinterpret_cast<const char*>(atmo->nsxinte),
            static_cast<std::streamsize>(intensity_count * sizeof(double)));
  out.close();
  if (!out) {
    std::remove(temporary_path.c_str());
    return;
  }

  // A concurrent process may have completed the same immutable cache first.
  // On Windows rename does not replace an existing destination; replacing it
  // with the same source-validated payload keeps the operation race-safe.
#ifdef _WIN32
  if (std::rename(temporary_path.c_str(), NSX_H_CACHE_PATH) != 0) {
    // Only this generated cache is replaced. Concurrent writers have the
    // same validated source identity and numerical payload.
    std::remove(NSX_H_CACHE_PATH);
    if (std::rename(temporary_path.c_str(), NSX_H_CACHE_PATH) != 0) {
      std::remove(temporary_path.c_str());
    }
  }
#else
  if (std::rename(temporary_path.c_str(), NSX_H_CACHE_PATH) != 0) {
    std::remove(temporary_path.c_str());
  }
#endif
}

} // namespace

bool NSXHTableExists() {
  FILE *Hspecttable = fopen(nsx_h_source_path(), "r");
  if (!Hspecttable) return false;
  fclose(Hspecttable);
  return true;
}



// Computes the correct theta index for NSX
int th_index_nsx(double cos_alpha, class Atmo* atmo){

  int n_mu;

  //std::cout << "th_index_nsx: cos_theta = " << cos_theta << std::endl;  
  //Find closest value of mu=cos(alpha) in the table
    n_mu = 1;
    while (cos_alpha < atmo->nsxangl[n_mu] && n_mu < atmo->Nmu){
    	n_mu += 1;
    }
   
    return (n_mu-1);
}



// Read in Wynn Ho's Hydrogen Table
void ReadNSXH(class Atmo* atmo){

	double loget,mu,logt,logg,logi;
	// loget = log_10(E/kT) where E = photon energy, kT = Thermal energy in same units
	// mu = cos(theta); theta = Angle measured from the normal to the surface
	// logt = log_10(T) where T is in Kelvin
	// logg = log_10(g) where g is acceleration due to gravity in cgs units
	// logi = log_10(I/T^3)
	
	int NlogTeff, Nlogg, NlogE, Nmu, Npts;

	NlogTeff = 35;
	atmo->NlogTeff = NlogTeff;
	Nlogg = 14; // updated 20221116
	atmo->Nlogg = Nlogg;
	NlogE = 166;
	atmo->NlogE = NlogE;
	Nmu = 67;
	atmo->Nmu = Nmu;
	
	atmo->nsxlogTeff = dvector(0,NlogTeff);
	for (int i=0;i<NlogTeff;i++){
	  atmo->nsxlogTeff[i] = 5.10 + i*0.05;	  
	}


	atmo->nsxlogg = dvector(0,Nlogg);
	for (int i=0;i<Nlogg;i++){
	  atmo->nsxlogg[i] = 13.7 + 0.1*i;
	}

	atmo->nsxloget = dvector(0,NlogE);
	for (int i=0;i<NlogE;i++){
	  atmo->nsxloget[i] = -1.30 + i*0.02;
	}


	Npts =  (NlogTeff*Nlogg*NlogE);

	atmo->Npts = Npts;

	//std::cout << "Npts = " << Npts << std::endl;

	atmo->nsxangl = dvector(0,Nmu);
	atmo->nsxinte = dvector(0,Npts*Nmu);

    NSXSourceIdentity source = {};
    const char* source_path = nsx_h_source_path();
    if (!source_identity(source_path, &source)) {
      throw std::runtime_error(std::string("ReadNSXH: cannot stat hydrogen atmosphere table: ") + source_path);
    }
    if (load_nsx_cache(atmo, source)) return;

    FILE *Hspecttable = fopen(source_path, "r");
    if (!Hspecttable) {
      throw std::runtime_error("ReadNSXH: cannot open hydrogen atmosphere table.");
    }

	//	std::cout << "Finished allocating memory " << std::endl;
	
    	for (int i = 0; i < Npts*Nmu; i++){
	  
	  const int fields = fscanf(Hspecttable,"%lf %lf %lf %lf %lf",
		 &loget, &mu,
		 &logi, &logt, &logg);
      if (fields != 5) {
        fclose(Hspecttable);
        throw std::runtime_error("ReadNSXH: malformed or truncated hydrogen atmosphere table.");
      }
 
	  atmo->nsxinte[i] = logi;
		    
	  if ( i < Nmu ) {
	    atmo->nsxangl[i] = mu;
	  }
	  
	  
	}
       
	fclose(Hspecttable);
	write_nsx_cache(atmo, source);
    	//std::cout << "finished reading Wynn Ho's NSX-H Atmosphere" << std::endl;
}



// NSX - Hydrogen Atmosphere Computed by Wynn Ho
// Given values of T, g, and emission angle, compute the function I(E), stored in spectrum
// Computation for "on-grid" values
// This "new" version takes into account that the energy is really the ratio: E/kT
void NSXHSpectrum(double cos_alpha, double lt, double lgrav, class Atmo* nsx, class Spectrum* spectrum){

  // lt = log_10(T) where T is in Kelvin
  // lgrav = log_10(g) where g is acceleration due to gravity in cgs units
  
  class Atmo atmo;
  atmo = (*nsx);

  int index;

  // Calculate the angle index closest to the given cos_alpha value
  int ii_mu = th_index_nsx(cos_alpha, &atmo);
  double dmu = cos_alpha - atmo.nsxangl[ii_mu];
  int kmmax(1);
  if (dmu!=0.0)
    kmmax=2;


  
  // Find the correct temperature index
  int i_lt = (lt-5.1)/0.05 ; //if we need to load 1st temperature, i_lt = 0. 
  double dlogt = lt - (i_lt*0.05 + 5.1);
  int ktmax(1);
  if (dlogt!=0.0)
    ktmax=2;

  // Find the correct gravity index
  int i_lgrav = (lgrav-13.7)/0.1;
  double dlogg = lgrav - (i_lgrav * 0.1 + 13.7);
  int kgmax(1);
  if (dlogg!=0.0)
    kgmax=2;
    

    double loget;

    double logi[3][3][3];
    double logim[3][3];
    double logig[3];
    double logit;

      
    for( int j(0); j<atmo.NlogE; j++){
      loget = atmo.nsxloget[j]; // log_10(E/T)

      for (int kt(0); kt<ktmax; kt++){ //loop through log(T) grid points
	for (int kg(0); kg<kgmax; kg++){ //loop through log(g) grid points
	  for (int km(0); km<kmmax; km++){ //loop through mu grid points
      
	      index = (i_lt+kt) * atmo.NlogE * atmo.Nmu * atmo.Nlogg
		+ (i_lgrav+kg) * atmo.NlogE * atmo.Nmu
		+ (j) * atmo.Nmu
		+ (ii_mu+km) ;
	    
	      logi[kt][kg][km] = atmo.nsxinte[index]; // log_10(I/T^3)	

	      /* if (j==50) // Checking values
		std::cout << "kt = " << kt
			<< " kg = " << kg
			<< " km = " << km
			<< " logi = " << logi[kt][kg][km]
			<< std::endl; */

	      
	  } // end km loop

	  if (kmmax!=1){ // interpolate to correct mu value
	    logim[kt][kg] = logi[kt][kg][0]
	      + (logi[kt][kg][0]-logi[kt][kg][1])/(atmo.nsxangl[ii_mu]-atmo.nsxangl[ii_mu+1])*dmu;
	  }
	  else logim[kt][kg] = logi[kt][kg][0];

	  /* if (j==50)
	    std::cout << " Interpolated Value = logim[kt][kg] = " << logim[kt][kg]
	    << std::endl;*/
	  
	} // end kg loop
	if (kgmax!=1){ // interpolate to correct logg value
	  logig[kt] = logim[kt][0]
	    + (logim[kt][0] - logim[kt][1])/(atmo.nsxlogg[i_lgrav]-atmo.nsxlogg[i_lgrav+1])*dlogg;
	}
	else logig[kt] = logim[kt][0];

	/*	if (j==50)
	  std::cout << " Interpolated Value = logig[kt] = " << logig[kt]
	  << std::endl;*/	
      } // end kt loop
      if (ktmax!=1){ // interpolate to correct logT value
	logit = logig[0]
	  + (logig[0]-logig[1])/(atmo.nsxlogTeff[i_lt]-atmo.nsxlogTeff[i_lt+1])*dlogt;
      }
      else
	logit = logig[0];

      /*if (j==50)
	std::cout << "               Interpolated Value = logit = " << logit << std::endl;*/

      // These should be the final correct interpolated values 
      spectrum->energy[j] = pow(10,loget+lt) * Units::K_BOLTZ / Units::EV *1e-3 ; // Photon energy in keV
      spectrum->intensity[j] =  pow(10.0,logit+lt*3.0); // Specific Intensity in units of:
      // (erg s^-1 cm^-2 Hz^-1 ster^-1)

    }
	 	
}

// Free allocated memory
void CleanNSXH(class Atmo* atmo, class Spectrum* spectrum){

  
       free_dvector(atmo->nsxlogTeff,0,atmo->NlogTeff);
       free_dvector(atmo->nsxlogg,0,atmo->Nlogg);
       free_dvector(atmo->nsxloget,0,atmo->NlogE);
       free_dvector(atmo->nsxangl,0,atmo->Nmu);
       free_dvector(atmo->nsxinte,0,atmo->Npts*atmo->Nmu);

       free_dvector(spectrum->energy,0,atmo->NlogE);
       free_dvector(spectrum->intensity,0,atmo->NlogE);
 
}

