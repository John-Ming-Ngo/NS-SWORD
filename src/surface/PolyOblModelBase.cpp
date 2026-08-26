// PolyOblModelBase.cpp
//
// (C) Coire Cadeau, 2007

// Source (C) Coire Cadeau 2007, all rights reserved.
//
// Permission is granted for private use only, and not
// distribution, either verbatim or of derivative works,
// in whole or in part.
//
// The code is not thoroughly tested or guaranteed for
// any particular use.


#include "../../include/PolyOblModelBase.h"
#include "../../include/Exception.h"
#include "../../include/MathFunctions.h"
#include <cmath>


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

static std::string coefficients_csv = "inputs/shape_functions/Model_Legendre_Coefficients.csv";

PolyOblModelBase::PolyOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits)
  : OblModelBase(req_nounits, mass_nounits, omega_nounits), zeta(zetaparam(mass_nounits, req_nounits)), eps(epsparam(omega_nounits, mass_nounits, req_nounits)) { }

PolyOblModelBase::PolyOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model )
  : OblModelBase(req_nounits, mass_nounits, omega_nounits), zeta(zetaparam(mass_nounits, req_nounits)), eps(epsparam(omega_nounits, mass_nounits, req_nounits)), model(model) {

    // Load in the csv
    std::ifstream file(coefficients_csv);
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + coefficients_csv);
    }

    // Set the coefficients accordingly
    bool found = false;
    std::getline(file, line); //Throw away the header
    while ((std::getline(file, line)) && (!found)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;

        // Split the line by commas
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }

        // Convert the first token (model number) to an integer
        int csv_model = std::stoi(tokens[0]);

        // Check if this is the correct model
        if (csv_model == model) {
            found = true;
            model_name = tokens[1];  // model name

            a_0_0 = std::stod(tokens[2]);
            a_0_1 = std::stod(tokens[3]);
            a_0_2 = std::stod(tokens[4]);
            a_2_0 = std::stod(tokens[5]);
            a_2_1 = std::stod(tokens[6]);
            a_2_2 = std::stod(tokens[7]);
            a_4_0 = std::stod(tokens[8]);
            a_4_1 = std::stod(tokens[9]);
            a_4_2 = std::stod(tokens[10]);
        }
    }

    if (!found) {
        throw std::runtime_error("Model " + std::to_string(model) + " not found in the CSV file.");
    }

    file.close();
  }

bool PolyOblModelBase::validModel(int model) {
    bool found = false;
    // Load in the csv
    std::ifstream file(coefficients_csv);
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + coefficients_csv);
    }
    // Set the coefficients accordingly
    std::getline(file, line); //Throw away the header

    while ((std::getline(file, line)) && (!found)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;

        // Split the line by commas
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }

        // Convert the first token (model number) to an integer
        int csv_model = std::stoi(tokens[0]);

        // Check if this is the correct model
        if (csv_model == model) {
            found = true;
        }
    }
    file.close();

    return found;
}

double PolyOblModelBase::a0() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  return double(a_0_0*eps + a_0_1*zeta*eps + a_0_2*eps*eps);
}

double PolyOblModelBase::a2() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  return double(a_2_0*eps + a_2_1*zeta*eps + a_2_2*eps*eps);
}

double PolyOblModelBase::a4() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());

  return double(a_4_0*eps + a_4_1*zeta*eps + a_4_2*eps*eps);
}

double PolyOblModelBase::R_at_costheta( const double& costheta )  {
  // Return R(theta) in "nounits".
  // note that the user supplies cos(theta) and not theta.
  return double(Req_nounits()*( 1.0 + a0()*P0(costheta) + a2()*P2(costheta) + a4()*P4(costheta) ) ); 
}

double PolyOblModelBase::Dtheta_R( const double& costheta )  {
  // Return dR(theta) / dtheta in "nounits".
  // note that the user supplies cos(theta) and not theta.
  if(costheta == 0.0 || fabs(costheta) == 1.0) return double(0.0);
  else
    // -sin(theta) * Dcostheta(R)
    return double( -sqrt(1.0-costheta*costheta)  // -sin(theta). Our theta runs from 0 to pi; so sin(theta) is always positive, so this is acceptable.
		   * Req_nounits()
		   *( a0()*Dmu_P0(costheta) + a2()*Dmu_P2(costheta) + a4()*Dmu_P4(costheta) )
		   );
}

double PolyOblModelBase::zetaparam( const double& Mass_nounits, const double& Req_nounits ) {
  return double( Mass_nounits / Req_nounits );
}

double PolyOblModelBase::epsparam( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits ) {
  return double( pow(Omega_nounits,2.0) * pow(Req_nounits,3.0) / Mass_nounits );
}

double PolyOblModelBase::get_zeta() const {
  return double(zeta);
}

double PolyOblModelBase::get_eps() const {
  return double(eps);
}
