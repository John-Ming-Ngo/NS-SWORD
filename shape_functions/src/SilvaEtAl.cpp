#include "SilvaEtAl.h"
#include "../../include/Exception.h"
#include "../../include/MathFunctions.h"
#include <cmath>


#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

const std::string coefficients_csv = "inputs/shape_functions/Silva_Eliptical_Coefficients.csv";

SilvaOblModel::SilvaOblModel(double req_nounits, double mass_nounits, double omega_nounits) 
    : OblModelBase(req_nounits, mass_nounits, omega_nounits), zeta(zetaparam(mass_nounits, req_nounits)), eps(epsparam(omega_nounits, mass_nounits, req_nounits)) {
    
    //std::cout <<
    //"Mass: " << mass_nounits <<
    //", Radius: " << req_nounits <<
    //", Omega: " << omega_nounits <<
    //", Spin Param: " << eps <<
    //std::endl;

    //exit(1);
    //Choose the model according to Silva's criteria with eps.
        // Helper lambda to trim whitespace from a string
    auto trim = [](std::string& str) {
        str.erase(0, str.find_first_not_of(" \t\n\r"));
        str.erase(str.find_last_not_of(" \t\n\r") + 1);
    };

    int model = (eps <= 0.2) ? 1 : 2;

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
            model_name = "Silva_" + tokens[1];  // model name

            e_c_0_0     = std::stod(tokens[2]);
            e_c_1_2_0   = std::stod(tokens[3]);
            e_c_1_0     = std::stod(tokens[4]);
            e_c_0_1     = std::stod(tokens[5]);
            e_c_1_1     = std::stod(tokens[6]);
            e_c_2_0     = std::stod(tokens[7]);
            e_c_0_2     = std::stod(tokens[8]);

            a2_c_0_0    = std::stod(tokens[9]);
            a2_c_1_2_0  = std::stod(tokens[10]);
            a2_c_1_0    = std::stod(tokens[11]);
            a2_c_0_1    = std::stod(tokens[12]);
            a2_c_1_1    = std::stod(tokens[13]);
            a2_c_2_0    = std::stod(tokens[14]);
            a2_c_0_2    = std::stod(tokens[15]);

            a4_c_0_0    = std::stod(tokens[16]);
            a4_c_1_2_0  = std::stod(tokens[17]);
            a4_c_1_0    = std::stod(tokens[18]);
            a4_c_0_1    = std::stod(tokens[19]);
            a4_c_1_1    = std::stod(tokens[20]);
            a4_c_2_0    = std::stod(tokens[21]);
            a4_c_0_2    = std::stod(tokens[22]);
        }
    }

    if (!found) {
        throw std::runtime_error("Model " + std::to_string(model) + " not found in the CSV file.");
    }

    file.close();
}

SilvaOblModel::SilvaOblModel(double req_nounits, double mass_nounits, double omega_nounits, int model) 
    : OblModelBase(req_nounits, mass_nounits, omega_nounits), zeta(zetaparam(mass_nounits, req_nounits)), eps(epsparam(omega_nounits, mass_nounits, req_nounits)) {
    // Load in the csv
    std::ifstream file(coefficients_csv);
    std::string line;

    // Helper lambda to trim whitespace from a string
    auto trim = [](std::string& str) {
        str.erase(0, str.find_first_not_of(" \t\n\r"));
        str.erase(str.find_last_not_of(" \t\n\r") + 1);
    };

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
            model_name = tokens[1];  // model name

            e_c_0_0     = std::stod(tokens[2]);
            e_c_1_2_0   = std::stod(tokens[3]);
            e_c_1_0     = std::stod(tokens[4]);
            e_c_0_1     = std::stod(tokens[5]);
            e_c_1_1     = std::stod(tokens[6]);
            e_c_2_0     = std::stod(tokens[7]);
            e_c_0_2     = std::stod(tokens[8]);

            a2_c_0_0    = std::stod(tokens[9]);
            a2_c_1_2_0  = std::stod(tokens[10]);
            a2_c_1_0    = std::stod(tokens[11]);
            a2_c_0_1    = std::stod(tokens[12]);
            a2_c_1_1    = std::stod(tokens[13]);
            a2_c_2_0    = std::stod(tokens[14]);
            a2_c_0_2    = std::stod(tokens[15]);

            a4_c_0_0    = std::stod(tokens[16]);
            a4_c_1_2_0  = std::stod(tokens[17]);
            a4_c_1_0    = std::stod(tokens[18]);
            a4_c_0_1    = std::stod(tokens[19]);
            a4_c_1_1    = std::stod(tokens[20]);
            a4_c_2_0    = std::stod(tokens[21]);
            a4_c_0_2    = std::stod(tokens[22]);
        }
    }

    if (!found) {
        throw std::runtime_error("Model " + std::to_string(model) + " not found in the CSV file.");
    }

    file.close();
}

SilvaOblModel::SilvaOblModel(double req_nounits, double mass_nounits, double omega_nounits, std::string model_name)
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
        if (tokens.size() < 23) {
            throw std::runtime_error("Malformed line in CSV file: " + line);
        }

        // Trim the model name token
        std::string csv_model_name = tokens[1];
        trim(csv_model_name);

        // Check if this is the correct model
        if (csv_model_name == model_name) {
            found = true;
            this->model_name = csv_model_name;  // Assign trimmed model name

            e_c_0_0     = std::stod(tokens[2]);
            e_c_1_2_0   = std::stod(tokens[3]);
            e_c_1_0     = std::stod(tokens[4]);
            e_c_0_1     = std::stod(tokens[5]);
            e_c_1_1     = std::stod(tokens[6]);
            e_c_2_0     = std::stod(tokens[7]);
            e_c_0_2     = std::stod(tokens[8]);

            a2_c_0_0    = std::stod(tokens[9]);
            a2_c_1_2_0  = std::stod(tokens[10]);
            a2_c_1_0    = std::stod(tokens[11]);
            a2_c_0_1    = std::stod(tokens[12]);
            a2_c_1_1    = std::stod(tokens[13]);
            a2_c_2_0    = std::stod(tokens[14]);
            a2_c_0_2    = std::stod(tokens[15]);

            a4_c_0_0    = std::stod(tokens[16]);
            a4_c_1_2_0  = std::stod(tokens[17]);
            a4_c_1_0    = std::stod(tokens[18]);
            a4_c_0_1    = std::stod(tokens[19]);
            a4_c_1_1    = std::stod(tokens[20]);
            a4_c_2_0    = std::stod(tokens[21]);
            a4_c_0_2    = std::stod(tokens[22]);
        }
    }

    if (!found) {
        throw std::runtime_error("Model " + model_name + " not found in the CSV file.");
    }

    file.close();
}

bool SilvaOblModel::validModel(int model) {
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

double SilvaOblModel::zetaparam( const double& Mass_nounits, const double& Req_nounits ) {
  return double( Mass_nounits / Req_nounits );
}

double SilvaOblModel::epsparam( const double& Omega_nounits, const double& Mass_nounits, const double& Req_nounits ) {
  return double( pow(Omega_nounits,2.0) * pow(Req_nounits,3.0) / Mass_nounits );
}

double SilvaOblModel::get_zeta() const {
  return double(zeta);
}

double SilvaOblModel::get_eps() const {
  return double(eps);
}

double SilvaOblModel::y(double c_0_0, double c_1_2_0, double c_1_0, double c_0_1, double c_1_1, double c_2_0, double c_0_2,
 double eps, double zeta) const {
  return c_0_0 + c_1_2_0 * pow(eps, 0.5) + c_1_0 * pow(eps, 1) + c_0_1 * zeta + c_1_1 * eps * zeta + c_2_0 * pow(eps, 2) + c_0_2 * pow(zeta, 2);
}

double SilvaOblModel::eccentricity() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  return double(y(e_c_0_0, e_c_1_2_0, e_c_1_0, e_c_0_1, e_c_1_1, e_c_2_0, e_c_0_2, eps, zeta));
}

double SilvaOblModel::a2() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());
  return double(y(a2_c_0_0, a2_c_1_2_0, a2_c_1_0, a2_c_0_1, a2_c_1_1, a2_c_2_0, a2_c_0_2, eps, zeta));
}

double SilvaOblModel::a4() const {
  double eps(this->get_eps());
  double zeta(this->get_zeta());

  return double(y(a4_c_0_0, a4_c_1_2_0, a4_c_1_0, a4_c_0_1, a4_c_1_1, a4_c_2_0, a4_c_0_2, eps, zeta));
}

double SilvaOblModel::g(const double& costheta) const {
  return (1 + a2() * pow(costheta, 2) + a4() * pow(costheta, 4) - (1 + a2() + a4()) * pow(costheta, 6));
}

double SilvaOblModel::dg_dcostheta(const double& costheta) const {
  return (a2() * 2 * pow(costheta, 1) 
  + a4() * 4 * pow(costheta, 3) 
  - (1 + a2() + a4()) * 6 * pow(costheta, 5));
}

double SilvaOblModel::R_at_costheta( const double& costheta ) const {
  // Return R(theta) in "nounits".
  // note that the user supplies cos(theta) and not theta.
  double req = Req_nounits();
  double radical = pow((1-pow(eccentricity(), 2))/(1-(pow(eccentricity(), 2))*g(costheta)), 0.5);
  //std::cout << "costheta: " << costheta << ", radical: " << radical 
  //<< ", Inside term: " << (1-pow(eccentricity(), 2))/(1-(pow(eccentricity(), 2))*g(costheta)) 
  //<< ", eccentricity:" << eccentricity() 
  //<< ", g(costheta)" << g(costheta) << std::endl;
  return req * radical; 
}

double SilvaOblModel::Dtheta_R( const double& costheta )  {
    // Return dR(theta) / dtheta in "nounits".
    // Note that the user supplies cos(theta) and not theta.
    double req = Req_nounits();
    double e2 = pow(eccentricity(), 2);
    double g_costheta = g(costheta);
    double radical = pow((1 - e2) / (1 - e2 * g_costheta), 0.5);

    // For q(mu) = sqrt((1-e^2)/(1-e^2 g(mu))),
    // dq/dmu = q e^2 g'(mu) / (2 (1-e^2 g(mu))).
    double dradical = radical * e2 * dg_dcostheta(costheta)
      / (2.0 * (1.0 - e2 * g_costheta));

    // Derivative of cos(theta) with respect to theta:
    double dcostheta_dtheta = -sqrt(1.0 - costheta * costheta); // -sin(theta). Our theta runs from 0 to pi; so sin(theta) is always positive, so this is acceptable.

    // Combine using the chain rule:
    return req * (dradical * dcostheta_dtheta);
}

extern "C" OblModelBase* createOblModelBase(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, int model) {
  return new SilvaOblModel(req_nounits, mass_nounits, omega_nounits, model); //Use the model assigned by the user.
  //return new SilvaOblModel(req_nounits, mass_nounits, omega_nounits); //We assign the model according to a program parameter ourselves. This follows from Silva's modelling criteria based on the dimensionless spin parameter cutoff.
}

extern "C" OblModelBase* createOblModelBaseStr(const double& req_nounits, const double& mass_nounits, const double& omega_nounits, std::string model_name) {
    return new SilvaOblModel(req_nounits, mass_nounits, omega_nounits, model_name);
}
