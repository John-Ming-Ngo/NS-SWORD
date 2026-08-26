// ElipticalMorsink.cpp
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


#include "ElipticalMorsink.h"
#include "../../include/Exception.h"
#include "../../include/MathFunctions.h"
#include "../../include/Units.h"
#include <cmath>


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

const std::string coefficients_csv = "inputs/shape_functions/Model_Legendre_Coefficients.csv";

ElipticalMorsink::ElipticalMorsink(const double& req_nounits, const double& mass_nounits, const double& omega_nounits)
  : OblModelBase(req_nounits, mass_nounits, omega_nounits), zeta(zetaparam(mass_nounits, req_nounits)), eps(epsparam(omega_nounits, mass_nounits, req_nounits)) { }

ElipticalMorsink::ElipticalMorsink(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model )
  : OblModelBase(req_nounits, mass_nounits, omega_nounits), zeta(zetaparam(mass_nounits, req_nounits)), eps(epsparam(omega_nounits, mass_nounits, req_nounits)), model(model) {
    // Helper lambda to trim whitespace from a string
    auto trim = [](std::string& str) {
      str.erase(0, str.find_first_not_of(" \t\n\r"));
      str.erase(str.find_last_not_of(" \t\n\r") + 1);
    };

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
            trim(tokens[1]);
            model_name = "Elliptical_Morsink_" + tokens[1];  // model name

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

ElipticalMorsink::ElipticalMorsink(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name)
  : OblModelBase(req_nounits, mass_nounits, omega_nounits), zeta(zetaparam(mass_nounits, req_nounits)), eps(epsparam(omega_nounits, mass_nounits, req_nounits)) {

    // Load in the csv
    std::ifstream file(coefficients_csv);
    std::string line;

    if (!file.is_open()) {
        throw std::runtime_error("Unable to open file: " + coefficients_csv);
    }

    // Helper lambda to trim whitespace from a string
    auto trim = [](std::string& str) {
        str.erase(0, str.find_first_not_of(" \t\n\r"));
        str.erase(str.find_last_not_of(" \t\n\r") + 1);
    };

    // Trim the input model name
    trim(model_name);

    // Set the coefficients accordingly
    bool found = false;
    std::getline(file, line); // Throw away the header
    while ((std::getline(file, line)) && (!found)) {
        std::stringstream ss(line);
        std::string token;
        std::vector<std::string> tokens;

        // Split the line by commas
        while (std::getline(ss, token, ',')) {
            tokens.push_back(token);
        }

        // Ensure there are enough tokens to process
        if (tokens.size() < 11) {
            throw std::runtime_error("Malformed line in CSV file: " + line);
        }

        // Trim the model name token
        std::string csv_model_name = tokens[1];
        trim(csv_model_name);

        // Check if this is the correct model
        if (csv_model_name == model_name) {
            found = true;
            this->model_name = "Elliptical_Morsink_" + csv_model_name;  // Assign trimmed model name

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
        throw std::runtime_error("Model " + model_name + " not found in the CSV file.");
    }

    file.close();
}

bool ElipticalMorsink::validModel(int model) {
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

double ElipticalMorsink::a0() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  return double(a_0_0*eps + a_0_1*zeta*eps + a_0_2*eps*eps);
}

double ElipticalMorsink::a2() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  return double(a_2_0*eps + a_2_1*zeta*eps + a_2_2*eps*eps);
}

double ElipticalMorsink::a4() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());

  return double(a_4_0*eps + a_4_1*zeta*eps + a_4_2*eps*eps);
}

double ElipticalMorsink::R_at_costheta( const double& costheta ) const {
  // Return R(theta) in "nounits".
  // note that the user supplies cos(theta) and not theta.
  // Error Diagnostic
  //std::cout 
  //  << "costheta: " << costheta 
  //  << ", theta: " << acos(costheta) 
  //  << ", radius calculated: " << double(Req_nounits()*( 1.0 + a0()*P0(costheta) + a2()*P2(costheta) + a4()*P4(costheta) ) )
  //  << ", original radius: " << Req_nounits()
  //  << ", ratio in unitless: " << double(Req_nounits()*( 1.0 + a0()*P0(costheta) + a2()*P2(costheta) + a4()*P4(costheta) ) )/Req_nounits()
  //  << std::endl;
  double theta = acos(costheta);
  double req = Req_nounits();
  double costheta_pole = cos(0);
  double rpolar = double(Req_nounits()*( 1.0 + a0()*P0(costheta_pole) + a2()*P2(costheta_pole) + a4()*P4(costheta_pole) ) );
  double ratio = rpolar/req;
  return req * (pow(sin(theta), 2) + ratio * pow(cos(theta), 2));
}

double ElipticalMorsink::Dtheta_R( const double& costheta ) {
  // Return dR(theta) / dtheta in "nounits".
  // note that the user supplies cos(theta) and not theta.
  double theta = acos(costheta);
  double req = Req_nounits();
  double costheta_pole = cos(0);
  double rpolar = double(Req_nounits()*( 1.0 + a0()*P0(costheta_pole) + a2()*P2(costheta_pole) + a4()*P4(costheta_pole) ) );
  double ratio = rpolar/req;
  return req * (2*sin(theta)*cos(theta) - ratio * 2*cos(theta)*sin(theta));
}

double ElipticalMorsink::zetaparam( const double& Mass_nounits, const double& Req_nounits ) {
  return double( Mass_nounits / Req_nounits );
}

double ElipticalMorsink::epsparam( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits ) {
  return double( pow(Omega_nounits,2.0) * pow(Req_nounits,3.0) / Mass_nounits );
}

double ElipticalMorsink::get_zeta() const {
  return double(zeta);
}

double ElipticalMorsink::get_eps() const {
  return double(eps);
}

extern "C" OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model) {
    return new ElipticalMorsink(req_nounits, mass_nounits, omega_nounits, model);
}

extern "C" OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name) {
    return new ElipticalMorsink(req_nounits, mass_nounits, omega_nounits, model_name);
}
