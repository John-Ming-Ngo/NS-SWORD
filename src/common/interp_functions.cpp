#include "../../include/interp_functions.h"
#include <cmath>

#include <map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <exception>

// Polynomial interpolation function using Lagrange method - taken from AI. I do not understand this yet. 
double polynomial_interpolate(const std::vector<double>& xp, const std::vector<double>& yp, double x, int order) {
    if (xp.size() != yp.size() || xp.empty() || order < 1) {
        throw std::invalid_argument("xp and yp must be of the same non-zero size, and order must be at least 1.");
    }
    if (order + 1 > xp.size()) {
        throw std::invalid_argument("Order is too high for the number of data points.");
    }

    // Extrapolation for out-of-bounds x
    if (x < xp.front()) {
        double slope = (yp[1] - yp[0]) / (xp[1] - xp[0]);
        return yp[0] + slope * (x - xp[0]);
    }
    if (x > xp.back()) {
        size_t n = xp.size();
        double slope = (yp[n - 1] - yp[n - 2]) / (xp[n - 1] - xp[n - 2]);
        return yp[n - 1] + slope * (x - xp[n - 1]);
    }

    // Find the closest point in xp to x
    auto closest = std::lower_bound(xp.begin(), xp.end(), x);
    size_t center = closest == xp.begin() ? 0 : closest - xp.begin() - 1;

    // Adjust indices to select `order + 1` points centered around `center`
    int half = order / 2;
    int start = std::max(0, static_cast<int>(center) - half);
    int end = std::min(static_cast<int>(xp.size()), start + order + 1);
    start = end - order - 1;

    // Calculate Lagrange interpolation polynomial
    double result = 0.0;
    for (int i = start; i < end; ++i) {
        double term = yp[i];
        for (int j = start; j < end; ++j) {
            if (i != j) {
                term *= (x - xp[j]) / (xp[i] - xp[j]);
            }
        }
        result += term;
    }

    return result;
}

// Linear interpolation function
double linear_interpolate(const std::vector<double>& xp, const std::vector<double>& yp, double x) {
    if (xp.size() != yp.size() || xp.empty()) {
        throw std::invalid_argument("xp and yp must be of the same non-zero size.");
    }

    // Extrapolation for out-of-bounds x
    if (x < xp.front()) {
        double slope = (yp[1] - yp[0]) / (xp[1] - xp[0]);
        return yp[0] + slope * (x - xp[0]);
    }
    if (x > xp.back()) {
        size_t n = xp.size();
        double slope = (yp[n - 1] - yp[n - 2]) / (xp[n - 1] - xp[n - 2]);
        return yp[n - 1] + slope * (x - xp[n - 1]);
    }

    // Find the interval [x0, x1] around x
    auto upper = std::lower_bound(xp.begin(), xp.end(), x);
    if (upper == xp.begin()) {
        return yp.front();
    }
    if (upper == xp.end()) {
        return yp.back();
    }

    // Indices of the two points for interpolation
    size_t i1 = upper - xp.begin();
    size_t i0 = i1 - 1;

    // Perform linear interpolation
    double x0 = xp[i0];
    double y0 = yp[i0];
    double x1 = xp[i1];
    double y1 = yp[i1];

    return y0 + (x - x0) * (y1 - y0) / (x1 - x0);
}
