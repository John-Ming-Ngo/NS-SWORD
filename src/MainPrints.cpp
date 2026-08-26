#include "MainPrints.h"

void PrintMainFlags(const MainFlags* flags, std::ostream& out)
{
    if (!flags) {
        std::cerr << "Error: Null pointer to MainFlags.\n";
    }

    out << "=== MainFlags ===\n";
    out << "shape_model_set: "               << (flags->shape_model_set             ? "true" : "false") << "\n";
    out << "temperature_model_set: "         << (flags->temperature_model_set       ? "true" : "false") << "\n";

    out << "mass_set: "                     << (flags->mass_set                    ? "true" : "false") << "\n";
    out << "req_set: "                      << (flags->req_set                     ? "true" : "false") << "\n";
    out << "Teq_set: "                      << (flags->Teq_set                     ? "true" : "false") << "\n";
    out << "omega_set: "                    << (flags->omega_set                   ? "true" : "false") << "\n";
    out << "inclination_set: "              << (flags->inclination_set             ? "true" : "false") << "\n";
    out << "distance_set: "                 << (flags->distance_set                ? "true" : "false") << "\n";
    out << "num_latitude_bins_set: "        << (flags->num_latitude_bins_set       ? "true" : "false") << "\n";
    out << "num_longitude_bins_set: "       << (flags->num_longitude_bins_set      ? "true" : "false") << "\n";

    out << "model_statistics_header_print: " << (flags->model_statistics_header_print ? "true" : "false") << "\n";

    out << "model_statistics_to_file: "      << (flags->model_statistics_to_file    ? "true" : "false") << "\n";
    out << "spectra_to_file: "              << (flags->spectra_to_file             ? "true" : "false") << "\n";
    out << "grid_to_file: "                 << (flags->grid_to_file                ? "true" : "false") << "\n";
}

bool PrintStarGridCoords(std::vector<std::vector<SpotData>> *StarGrid, std::ostream &out) {
  int ntheta = StarGrid->size();  // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      out << (*StarGrid)[i][j].latitude_bins_coord << ","
          << (*StarGrid)[i][j].longitude_bins_coord << ","
          << (*StarGrid)[i][j].theta << ","
          << (*StarGrid)[i][j].phi << ","
          << std::endl;
    }
  }
  return true;
}

bool PrintSpotData(std::vector<std::vector<SpotData>> *StarGrid, 
                   std::vector<std::string> print_header, 
                   std::string separator, 
                   std::ostream &out) {
  int ntheta = StarGrid->size();  // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (std::string value:print_header) {
    out << value << separator;
  }
  out << std::endl;

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &Spot = (*StarGrid)[i][j];
      for (std::string column : print_header) {
        out << Spot.DerivedVals[column] << separator;
      }
      out << std::endl;
    }
  }
  return true;
}

bool PrintOutput(std::map<std::string, std::string> outputs, 
                 std::vector<std::string> print_header, 
                 std::string separator, 
                 std::ostream &out) {
  bool found_all_columns = true;
  for (const std::string &key : print_header) {
    if (outputs.find(key) == outputs.end()) {
      found_all_columns = false;
      out << "" << separator;
    } else {
      out << outputs[key] << separator;
    }
  }
  out << std::endl;
  return found_all_columns;
}
