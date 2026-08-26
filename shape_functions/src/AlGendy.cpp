#include "AlGendy.h"
#include "../../include/Exception.h"
#include "../../include/MathFunctions.h"
#include <cmath>


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

const std::string coefficients_csv = "inputs/shape_functions/AlGendy_Eliptical_Coefficients.csv";

AlGendyOblModel::AlGendyOblModel(double req_nounits, double mass_nounits, double omega_nounits, int model) 
    : OblModelBase(req_nounits, mass_nounits, omega_nounits), zeta(zetaparam(mass_nounits, req_nounits)), eps(epsparam(omega_nounits, mass_nounits, req_nounits)) {
    
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
    std::getline(file, line); // Throw away the header
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
            model_name = "AlGendy_" + tokens[1];  // model name
            c_2_0 = std::stod(tokens[2]);
            c_2_1 = std::stod(tokens[3]);
        }
    }

    if (!found) {
        throw std::runtime_error("Model " + std::to_string(model) + " not found in the CSV file.");
    }

    file.close();
}

AlGendyOblModel::AlGendyOblModel(double req_nounits, double mass_nounits, double omega_nounits, std::string model_name) 
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
        if (tokens.size() < 4) {
            throw std::runtime_error("Malformed line in CSV file: " + line);
        }

        // Trim the model name token
        std::string csv_model_name = tokens[1];
        trim(csv_model_name);

        // Check if this is the correct model
        if (csv_model_name == model_name) {
            found = true;
            this->model_name = "AlGendy_" + csv_model_name;  // Assign trimmed model name
            c_2_0 = std::stod(tokens[2]);
            c_2_1 = std::stod(tokens[3]);
        }
    }

    if (!found) {
        throw std::runtime_error("Model " + model_name + " not found in the CSV file.");
    }

    file.close();
}

bool AlGendyOblModel::validModel(int model) {
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

double AlGendyOblModel::zetaparam( const double& Mass_nounits, const double& Req_nounits ) {
  return double( Mass_nounits / Req_nounits );
}

double AlGendyOblModel::epsparam( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits ) {
  return double( pow(Omega_nounits,2.0) * pow(Req_nounits,3.0) / Mass_nounits );
}

double AlGendyOblModel::get_zeta() const {
  return double(zeta);
}

double AlGendyOblModel::get_eps() const {
  return double(eps);
}

double AlGendyOblModel::a2() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  return double(c_2_0 * eps + c_2_1 * eps*zeta);
}

double AlGendyOblModel::R_at_costheta( const double& costheta ) const {
  // Return R(theta) in "nounits".
  // note that the user supplies cos(theta) and not theta.
  double req = Req_nounits();
  return req * (1 + a2()*pow(costheta, 2)); 
}

double AlGendyOblModel::Dtheta_R( const double& costheta )  {
    // Return dR(theta) / dtheta in "nounits".
    // Note that the user supplies cos(theta) and not theta.
    double req = Req_nounits();

    // Derivative of cos(theta) with respect to theta:
    double dcostheta_dtheta = -sqrt(1.0 - costheta * costheta); // -sin(theta). Our theta runs from 0 to pi; so sin(theta) is always positive, so this is acceptable.

    // Combine using the chain rule:
    return req * (2 * a2()*pow(costheta, 1) * dcostheta_dtheta);
}

extern "C" OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model) {
    return new AlGendyOblModel(req_nounits, mass_nounits, omega_nounits, model);
}

extern "C" OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name) {
    return new AlGendyOblModel(req_nounits, mass_nounits, omega_nounits, model_name);
}
