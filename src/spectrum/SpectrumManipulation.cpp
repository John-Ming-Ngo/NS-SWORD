#include <cmath>
#include <map>
#include <unistd.h>
#include "../../include/Hydrogen.h"
#include "../../include/BlackBody.h"
#include "../../include/Units.h"
#include "../../include/Struct.h"
#include "../../include/interp.h"
#include "../../include/nrutil.h"
#include <stdio.h>
#include <iostream>
#include <mutex>

using namespace std;

void SpectrumToMap(Spectrum *spectrum, std::map<double, double> *spectrumMap);

void PrintSpectrum(Spectrum spectrum) {
    for (int i = 0; i < spectrum.Npts; i++) {
        std::cout << spectrum.energy[i] << "," << spectrum.intensity[i] << std::endl;
    }
}

namespace {

const int NSX_ENERGY_POINTS = 166;
const double NSX_LOG_ENERGY_MIN = -1.30;
const double NSX_LOG_ENERGY_STEP = 0.02;

void InitializeEnergyGrid(Spectrum *spectrum, double logT) {
  // Allocate memory for spectrum
  spectrum->energy = dvector(0, NSX_ENERGY_POINTS);
  spectrum->intensity = dvector(0, NSX_ENERGY_POINTS);
  spectrum->Npts = NSX_ENERGY_POINTS;

  for (int i = 0; i < NSX_ENERGY_POINTS; ++i) {
    const double loget = NSX_LOG_ENERGY_MIN + i * NSX_LOG_ENERGY_STEP;
    spectrum->energy[i] = pow(10.0, loget + logT) * Units::K_BOLTZ / Units::EV * 1e-3;
    spectrum->intensity[i] = 0.0;
  }
}

void ReleaseSpectrum(Spectrum *spectrum) {
  if (!spectrum || spectrum->Npts <= 0) return;
  free_dvector(spectrum->energy, 0, NSX_ENERGY_POINTS);
  free_dvector(spectrum->intensity, 0, NSX_ENERGY_POINTS);
  spectrum->energy = NULL;
  spectrum->intensity = NULL;
  spectrum->Npts = 0;
}

typedef std::map<double, std::map<double, double> > SpectrumCache;

std::map<double, double> cached_empty_grid(double TkeV) {
  static SpectrumCache cache;
  static std::mutex cache_mutex;
  std::lock_guard<std::mutex> lock(cache_mutex);
  SpectrumCache::const_iterator found = cache.find(TkeV);
  if (found != cache.end()) return found->second;

  const double logT = log10(TkeV * 1000 * Units::EV_TO_K);
  Spectrum spectrum = {};
  InitializeEnergyGrid(&spectrum, logT);
  std::map<double, double> result;
  SpectrumToMap(&spectrum, &result);
  ReleaseSpectrum(&spectrum);
  cache[TkeV] = result;
  return result;
}

std::map<double, double> cached_blackbody(double TkeV) {
  static SpectrumCache cache;
  static std::mutex cache_mutex;
  std::lock_guard<std::mutex> lock(cache_mutex);
  SpectrumCache::const_iterator found = cache.find(TkeV);
  if (found != cache.end()) return found->second;

  const double logT = log10(TkeV * 1000 * Units::EV_TO_K);
  Spectrum spectrum = {};
  InitializeEnergyGrid(&spectrum, logT);
  BBSpectrum(logT, &spectrum);
  std::map<double, double> result;
  SpectrumToMap(&spectrum, &result);
  ReleaseSpectrum(&spectrum);
  cache[TkeV] = result;
  return result;
}

} // namespace

void InitializeSpectrum(Spectrum *spectrum, Atmo* atmo, double logT, double lgrav) {
  (void)atmo;
  (void)lgrav;
  InitializeEnergyGrid(spectrum, logT);
}

void SpectrumToMap(Spectrum *spectrum, std::map<double, double> *spectrumMap) {
    for (int i = 0; i < spectrum->Npts; i++) {
        // Check if the energy value is already in the map
        if (spectrumMap->find(spectrum->energy[i]) == spectrumMap->end()) {
            // If not, add the energy and intensity pair to the map
            (*spectrumMap)[spectrum->energy[i]] = spectrum->intensity[i];
        }
        //else {
        //    (*spectrumMap)[spectrum->energy[i]] += spectrum->intensity[i];
        //}
        // If the key already exists, do nothing (i.e., do not override)
    }
}

Atmo* obtainAtmo() {
    static class Atmo atmo;
    // C++11 guarantees one-time, thread-safe initialization of local statics.
    static const bool isAtmoLoaded = (ReadNSXH(&atmo), true);
    (void)isAtmoLoaded;
    return &atmo;
}

void GetHydrogenSpectrum(Spectrum *spectrum, double cos_alpha, double logT, double lgrav) {
    Atmo* atmo = obtainAtmo();
    InitializeEnergyGrid(spectrum, logT);
    NSXHSpectrum(cos_alpha, logT, lgrav, atmo, spectrum);
}

void GetBlackbodySpectrum(Spectrum *spectrum, double cos_alpha, double logT, double lgrav) {
    (void)cos_alpha;
    (void)lgrav;
    InitializeEnergyGrid(spectrum, logT);
    BBSpectrum(logT, spectrum);
}

double CalculateLGrav(double m_msun, double req_km) {
    double mass_over_req = m_msun/(req_km) * Units::GMC2 * 1;
    double req_cm = req_km * 1e5; // Convert to cm!

    double delta = 1.0 / sqrt(1.0 - 2.0*mass_over_req);
    double lgrav = log10(delta * mass_over_req * pow(Units::C,2)/req_cm);
    return lgrav;
}

void GetHSpectrumMap(std::map<double, double> *HMap, double cos_alpha, double TkeV, double m_msun, double req_km) {
    double lgrav = CalculateLGrav(m_msun, req_km);
    double logT = log10(TkeV * 1000 * Units::EV_TO_K);

    Spectrum HSpectrum;
    //std::cout << cos_alpha << "," << logT << "," << lgrav << std::endl;
    GetHydrogenSpectrum(&HSpectrum, cos_alpha, logT, lgrav);
    SpectrumToMap(&HSpectrum, HMap); //Transfers the info to the map.
    ReleaseSpectrum(&HSpectrum);
}

void GetBBSpectrumMap(std::map<double, double> *BBMap, double cos_alpha, double TkeV, double m_msun, double req_km) {
    (void)cos_alpha;
    (void)m_msun;
    (void)req_km;
    *BBMap = cached_blackbody(TkeV);
}

void InitializeSpectrumMap(std::map<double, double> *spectrumMap, double TkeV, double m_msun, double req_km) {
    (void)m_msun;
    (void)req_km;
    *spectrumMap = cached_empty_grid(TkeV);
}
