#include <cmath>
#include <map>
#include <vector>
#include <algorithm>
#include "../../include/Units.h"
#include "../../include/Struct.h"
#include "../../include/interp_functions.h" // New header for polynomial_interpolate
#include <iostream>
#include <stdio.h>

void spectrum_interpolate_add(std::map<double, double>* spectrum, std::map<double, double>* dspectrum, int window_size) {
    // Aggregated cell spectra already live on the run spectrum's exact energy
    // grid.  In that common case Lagrange interpolation at each identical node
    // returns the source value, so add it directly without allocating hundreds
    // of temporary interpolation vectors.
    if (spectrum->size() == dspectrum->size()) {
        bool identical_grid = true;
        auto destination = spectrum->cbegin();
        auto source = dspectrum->cbegin();
        for (; destination != spectrum->cend(); ++destination, ++source) {
            if (destination->first != source->first) {
                identical_grid = false;
                break;
            }
        }
        if (identical_grid) {
            auto destination_write = spectrum->begin();
            source = dspectrum->cbegin();
            for (; destination_write != spectrum->end(); ++destination_write, ++source) {
                destination_write->second += source->second;
            }
            return;
        }
    }

    // The source map is immutable during interpolation.  Convert it to
    // contiguous storage once; the previous implementation rebuilt the same
    // temporary vectors for every destination energy bin.
    std::vector<double> source_x;
    std::vector<double> source_y;
    source_x.reserve(dspectrum->size());
    source_y.reserve(dspectrum->size());
    for (const auto& point : *dspectrum) {
        source_x.push_back(point.first);
        source_y.push_back(point.second);
    }
    std::vector<double> xp;
    std::vector<double> yp;
    xp.reserve(static_cast<size_t>(2 * window_size + 1));
    yp.reserve(static_cast<size_t>(2 * window_size + 1));

    // Iterate over each point in the original spectrum
    for (auto& [x, y] : *spectrum) {
        const std::vector<double>::const_iterator lower =
            std::lower_bound(source_x.cbegin(), source_x.cend(), x);
        const size_t lower_index = static_cast<size_t>(lower - source_x.cbegin());
        const size_t start = lower_index > static_cast<size_t>(window_size)
            ? lower_index - static_cast<size_t>(window_size) : 0;
        int num_points = static_cast<int>(std::min(
            source_x.size() - start, static_cast<size_t>(2 * window_size + 1)));
        if (num_points % 2 == 0) { // Ensure odd number of points
            num_points--;
        }
        xp.resize(static_cast<size_t>(num_points));
        yp.resize(static_cast<size_t>(num_points));
        for (int i = 0; i < num_points; ++i) {
            xp[static_cast<size_t>(i)] = source_x[start + static_cast<size_t>(i)];
            yp[static_cast<size_t>(i)] = source_y[start + static_cast<size_t>(i)];
        }

        // Perform polynomial interpolation
        int order = std::min(3, num_points - 1); // Define polynomial order, e.g., cubic (adjustable)
        double interpolated_value = polynomial_interpolate(xp, yp, x, order);

        // Add the interpolated value to the original spectrum bin
        y += interpolated_value;
    }
}

double spectrum_integrate(const std::map<double, double>& spectrum, int order) {
    if (spectrum.size() < 2) {
        return 0.0; // Not enough points to integrate
    }

    double integral = 0.0;
    auto it = spectrum.begin();
    auto prev_it = it++;

    while (it != spectrum.end()) {
        double x0 = prev_it->first;
        double x1 = it->first;
        double y0 = prev_it->second;
        double y1 = it->second;

        double width = x1 - x0;
        double area = 0.0;

        if (order == 1) {
            // Trapezoidal rule
            area = 0.5 * (y0 + y1) * width;
        } else {
            // Simpson's rule or higher-order polynomial
            int n = std::min(order, static_cast<int>(spectrum.size()));
            std::vector<double> xp(n);
            std::vector<double> yp(n);

            xp[0] = x0;
            yp[0] = y0;
            xp[1] = x1;
            yp[1] = y1;

            // Fill xp and yp vectors as needed with points around the current segment
            auto it_prev = prev_it;
            auto it_next = it;

            for (int i = 2; i < n; ++i) {
                if (++it_next != spectrum.end()) {
                    xp[i] = it_next->first;
                    yp[i] = it_next->second;
                } else if (it_prev != spectrum.begin()) {
                    --it_prev;
                    xp[i] = it_prev->first;
                    yp[i] = it_prev->second;
                } else {
                    n = i; // fewer points than order
                    break;
                }
            }

            // Polynomial interpolation for midpoint
            double xb = (x0 + x1) / 2.0;
            double yb = polynomial_interpolate(xp, yp, xb, n - 1);

            // Simpson's rule approximation
            area = (width / 6.0) * (y0 + 4 * yb + y1);
        }

        integral += area;
        prev_it = it++;
    }
    return integral;
}

//===== Other Utility Functions =====//

double PlanckFunction(double T, double E) {
    double h_squared_c_squared = pow(Units::H_PLANCK, 2) * pow(Units::C, 2);
    double E_cubed = pow(E * 1000, 3) * pow(Units::EV, 3);
    return (2.0 * E_cubed / h_squared_c_squared) / (exp(E / T) - 1);
}

void spectrum_redshift(std::map<double, double>* spectrum, double redshift) {
    std::map<double, double> redshifted_spectrum;
    for (const auto& kv : *spectrum) {
        double new_key = kv.first * redshift;
        double new_value = kv.second * (redshift * redshift * redshift);
        redshifted_spectrum[new_key] = new_value;
    }
    spectrum->swap(redshifted_spectrum);
}

void spectrum_multiply(std::map<double, double>* spectrum, double factor) {
    for (auto& pair : *spectrum) {
        pair.second *= factor;
    }
}

void ReportSpectrumMap(std::map<double, double> map) {
    for(auto it = map.cbegin(); it != map.cend(); ++it) {
        std::cout << it->first << "," << it->second << std::endl;
    }
}
