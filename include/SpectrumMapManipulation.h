#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include "Units.h"
#include "Struct.h"
#include "interp.h"
#include <stdio.h>

void spectrum_interpolate_add(std::map<double, double>* spectrum, std::map<double, double>* dspectrum, int window_size);

double spectrum_integrate(const std::map<double, double>& spectrum, int order = 3);

double spectrum_redshift(std::map<double, double>* spectrum, double redshift);

void spectrum_multiply(std::map<double, double>* spectrum, double factor);

void ReportSpectrumMap(std::map<double, double> map);
