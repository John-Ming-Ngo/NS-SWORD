/* ---------- Papigkiotis.cpp ----------------------- */
#include "Papigkiotis.h"

#include <cstdio>          // popen, pclose
#include <cstdlib>         // std::strtod
#include <map>
#include <vector>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <iomanip>

#include "../../include/Exception.h"
#include "../../include/MathFunctions.h"
#include "../../include/Units.h"
#include "../../include/interp.h"
#include "../../include/interp_functions.h"

/*
 * Our PVLS implementation
 * Refactored in 2026 using ChatGPT. 
 * Note: Feature scaling is very important.
 */
#ifndef PAPIGKIOTIS_SURFACE_SCRIPT
#define PAPIGKIOTIS_SURFACE_SCRIPT "shape_functions/dependencies/Papigkiotis/ns_radius.py"
#endif

#ifndef PAPIGKIOTIS_DERIVATIVE_SCRIPT
#define PAPIGKIOTIS_DERIVATIVE_SCRIPT "shape_functions/dependencies/Papigkiotis/ns_log_derivative.py"
#endif

#ifndef PAPIGKIOTIS_MODEL_NAME
#define PAPIGKIOTIS_MODEL_NAME "PVLS"
#endif

#ifndef PAPIGKIOTIS_PYTHON
#ifdef _WIN32
#define PAPIGKIOTIS_PYTHON "python"
#else
#define PAPIGKIOTIS_PYTHON "python3"
#endif
#endif

#include "PVLSCache.h"

double Papigkiotis::C( const double& Mass_nounits, const double& Req_nounits ) {
  return double( Mass_nounits / Req_nounits );
}
double Papigkiotis::sigma( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits ) {
  return double( pow(Omega_nounits,2.0) * pow(Req_nounits,3.0) / Mass_nounits );
}

static std::map<double, double> call_python_derivative(const double Req, const double M, const double Omega)
{
    const double C     = Papigkiotis::C(M, Req);
    const double sigma = Papigkiotis::sigma(Omega, M, Req);

    std::ostringstream cmd;
    cmd << std::setprecision(17)
        << PAPIGKIOTIS_PYTHON " " PAPIGKIOTIS_DERIVATIVE_SCRIPT " "
        << Req << ' ' << C << ' ' << sigma;
    
    std::string full_cmd = cmd.str();
    std::ostringstream output_log;

    FILE* pipe = popen(cmd.str().c_str(), "r");
    if (!pipe)
        throw std::runtime_error("Papigkiotis: popen() failed for derivative.");

    std::map<double,double> muD;
    char buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe))
    {
        output_log << buffer;

        char* endptr = nullptr;
        double mu    = std::strtod(buffer, &endptr);

        // skip over any spaces/tabs before the comma
        while (endptr && (*endptr == ' ' || *endptr == '\t')) {
            ++endptr;
        }
        if (endptr && *endptr == ',') {
            // skip the comma and any following spaces
            ++endptr;
            while (*endptr == ' ') ++endptr;
            double dlogR = std::strtod(endptr, nullptr);
            muD.emplace(mu, dlogR);
        }
    }
    pclose(pipe);
    if (muD.empty()) {
        // dump all the context to stderr
        std::string all_output = output_log.str();  
        std::cerr
            << "Papigkiotis::call_python_derivative FAILED\n"
            << "  Req    = " << Req << "\n"
            << "  C      = " << C << "\n"
            << "  sigma  = " << sigma << "\n"
            << "  command= " << full_cmd << "\n"
            << "  output = \n" << all_output << "\n";
        throw std::runtime_error("Papigkiotis: no derivative data parsed.");
    }

    if (!muD.empty() && muD.begin()->first > -1.0) {
        for (auto& kv : muD) {
            double mu = kv.first;
            double dlogR = kv.second;
            if (mu > 0.0) {
                muD.emplace(-mu, -dlogR);
            }
        }
    }

    return muD;
}

static std::map<double,double> call_python_surface(const double Req, const double M, const double Omega)
{
    /* Build the C & σ inputs exactly as the Python expects */
    const double C      = Papigkiotis::C(M, Req);
    const double sigma  = Papigkiotis::sigma(Omega, M, Req);

    /* Command : python <script> R_eq C sigma */
    std::ostringstream cmd;
    cmd << std::setprecision(17)
        << PAPIGKIOTIS_PYTHON " " PAPIGKIOTIS_SURFACE_SCRIPT " "
        << Req << ' ' << C << ' ' << sigma;

    std::string full_cmd = cmd.str();
    std::ostringstream output_log;

    /* Launch and read back – POSIX popen() for brevity         */
    FILE* pipe = popen(cmd.str().c_str(), "r");
    if (!pipe)
        throw std::runtime_error("Papigkiotis: popen() failed.");

    std::map<double,double> muR;
    char  buffer[256];
    while (fgets(buffer, sizeof(buffer), pipe))
    {
        // Expected line: "<mu> , <R>\n"
        output_log << buffer;           // accumulate raw output for debug
        char* endptr = nullptr;
        double mu = std::strtod(buffer, &endptr);

        // skip over any spaces or tabs before the comma
        while (endptr && (*endptr == ' ' || *endptr == '\t')) {
            ++endptr;
        }
        if (endptr && *endptr == ',') {
            // advance past comma and any following space
            ++endptr;
            while (*endptr == ' ') ++endptr;
            double R = std::strtod(endptr, nullptr);
            muR.emplace(mu, R);
        }
    }
    pclose(pipe);
    if (muR.empty()) {
        // dump all the context to stderr
        std::string all_output = output_log.str();  
        std::cerr
            << "Papigkiotis::call_python_surface FAILED\n"
            << "  Req     = " << Req   << "\n"
            << "  C       = " << C     << "\n"
            << "  sigma   = " << sigma << "\n"
            << "  command = " << full_cmd << "\n"
            << "  output  =\n" << output_log.str() << "\n";
        throw std::runtime_error("Papigkiotis: no surface data parsed.");
    }

    if (!muR.empty() && muR.begin()->first > -1.0) {
        for (auto& kv : muR) {
            double mu = kv.first;
            double R = kv.second;
            if (mu > 0.0) {
                muR.emplace(-mu, R);
            }
        }
    }

    return muR;
}

/* ---------- constructor ------------------------------------- */
Papigkiotis::Papigkiotis(double req_nounits,
                         double mass_nounits,
                         double omega_nounits,
                         int model)
: OblModelBase(req_nounits, mass_nounits, omega_nounits)
{
    model_name = PAPIGKIOTIS_MODEL_NAME;

    const double C      = Papigkiotis::C(mass_nounits, req_nounits);
    const double sigma  = Papigkiotis::sigma(omega_nounits, mass_nounits, req_nounits);
    PVLSCache::require_trained_domain(C, sigma);
    //std::cout << C << "," << sigma << "," << "req_nounits" << std::endl;
    /* Cache the shape (μ, R) table for quick interpolation     */
    if (!PVLSCache::load(req_nounits, C, sigma,
                         &mu_radii_map, &mu_radii_derivative_map)) {
        mu_radii_map = call_python_surface(req_nounits,
                                           mass_nounits,
                                           omega_nounits);
        mu_radii_derivative_map = call_python_derivative(req_nounits,
                                                         mass_nounits,
                                                         omega_nounits);
    }
}

/* The overload taking a string forwards to the same helper */
Papigkiotis::Papigkiotis(double req_nounits,
                         double mass_nounits,
                         double omega_nounits,
                         std::string model_name)
: Papigkiotis(req_nounits, mass_nounits, omega_nounits, 0) {}

/* ---------- R(μ) via polynomial interpolation --------------- */


double Papigkiotis::R_at_costheta(const double& costheta) const
{
    if (mu_radii_map.empty())
        throw std::runtime_error("Papigkiotis::R_at_costheta - data not loaded.");

    /* Unpack map into two vectors for the helper routine       */
    std::vector<double> xp, yp;
    for (const auto& kv : mu_radii_map) {
        xp.push_back(kv.first);
        yp.push_back(kv.second);
    }
    const int order = 3;   // cubic polynomial (stable + cheap)
    return polynomial_interpolate(xp, yp, costheta, order);
}
/**
double Papigkiotis::Dtheta_R(const double& costheta) 
{
    // Derivative dR/dθ via numerical differentiation of R(μ)=R(cosθ)
    if (costheta == 0.0 || fabs(costheta) == 1.0) {
        return 0.0;
    }

    // 1) build vectors of (μ, R) from the cached map
    std::vector<double> xp, yp;
    xp.reserve(mu_radii_map.size());
    yp.reserve(mu_radii_map.size());
    for (const auto& kv : mu_radii_map) {
        xp.push_back(kv.first);
        yp.push_back(kv.second);
    }

    // 2) choose finite‐difference step in μ = cosθ
    const double epsilon = 1e-5;
    const int    order   = 3;  // same cubic interpolator

    // 3) interpolate R at μ, μ+ε, μ–ε
    double R0 = polynomial_interpolate(xp, yp, costheta,     order);
    double Rp = polynomial_interpolate(xp, yp, costheta+epsilon, order);
    double Rm = polynomial_interpolate(xp, yp, costheta-epsilon, order);

    // 4) dR/dμ ≈ (R(μ+ε) − R(μ−ε)) / (2ε)
    double dR_dmu = (Rp - Rm) / (2.0 * epsilon);

    // 5) convert to dR/dθ = (dR/dμ) * (dμ/dθ) = −sinθ * dR/dμ
    double sin_theta = std::sqrt(1.0 - costheta * costheta);
    return -sin_theta * dR_dmu;
}


 */
double Papigkiotis::Dtheta_R( const double& costheta )  {
    if (costheta == 0.0 || fabs(costheta) == 1.0) {
        return 0.0;
    }

    std::vector<double> xp, yp;
    xp.reserve(mu_radii_derivative_map.size());
    yp.reserve(mu_radii_derivative_map.size());
    for (const auto& kv : mu_radii_derivative_map) {
        xp.push_back(kv.first);   
        yp.push_back(kv.second);  
    }
    const int order = 3;
    double dLogR_dtheta = polynomial_interpolate(xp, yp, costheta, order);
    double dR_dtheta = dLogR_dtheta * R_at_costheta(costheta);
    //std::cout << costheta << "," << dR_dtheta << std::endl;
    //double sin_theta = std::sqrt(1.0 - costheta * costheta);
    return dR_dtheta;
}


extern "C" OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model) {
    return new Papigkiotis(req_nounits, mass_nounits, omega_nounits, model);
}

extern "C" OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name) {
    return new Papigkiotis(req_nounits, mass_nounits, omega_nounits, model_name);
}
