#include <cmath>

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <exception>

double polynomial_interpolate(const std::vector<double>& xp, const std::vector<double>& yp, double x, int order);
double linear_interpolate(const std::vector<double>& xp, const std::vector<double>& yp, double x);