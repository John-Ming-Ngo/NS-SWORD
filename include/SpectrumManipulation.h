#include <map>
#include <cmath>

#include "Hydrogen.h"
#include "BlackBody.h"
#include "Units.h"
#include "Struct.h"
#include <stdio.h>

using namespace std;

void InitializeSpectrum(Spectrum *spectrum, Atmo* atmo, double logT, double lgrav);

void SpectrumToMap(Spectrum *spectrum, std::map<double, double> *spectrumMap);

Atmo* obtainAtmo();

void GetHydrogenSpectrum(Spectrum *spectrum, double cos_alpha, double logT, double lgrav);

void GetBlackbodySpectrum(Spectrum *spectrum, double cos_alpha, double logT, double lgrav);

double CalculateLGrav(double m_msun, double req_km);

void GetHSpectrumMap(std::map<double, double> *HMap, double cos_alpha, double TkeV, double m_msun, double req_km);

void GetBBSpectrumMap(std::map<double, double> *BBMap, double cos_alpha, double TkeV, double m_msun, double req_km);

void InitializeSpectrumMap(std::map<double, double> *spectrumMap, double TkeV, double m_msun, double req_km);
