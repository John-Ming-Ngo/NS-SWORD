#ifndef MAINSTRUCTS_H
#define MAINSTRUCTS_H

#include <algorithm>
#include <iostream>
#include <fstream>
#include <cmath>
#include <exception>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>

#if defined(_WIN32) || defined(_WIN64)
#include <windows.h>
#else
#include <dlfcn.h>
#endif

#include <stdexcept>

#include "../include/OblDeflectionTOA.h"
#include "../include/Chi.h"
#include "../include/OblModelBase.h"
#include "../include/Units.h"
#include "../include/Exception.h"
#include "../include/Struct.h"

#include "../include/Hydrogen.h"
#include "../include/BlackBody.h"
#include "../include/SpectrumManipulation.h"
#include "../include/SpectrumMapManipulation.h"
#include <cstring> // For strcmp

//SpotData consists of everything one would want to know about a patch of the surface of the neutron star.
//This will be the primary data structure data transforms are performed on, once the structure is made.
struct SpotData
{
  //Global spot data - unclear if I should have anything for this. Ok, no, I shouldn't.

  //Local spot data
  double latitude_bins_coord; // Theta location on the theta-phi grid
  double theta; // Radians
  double longitude_bins_coord; // Phi location on the theta-phi grid
  double phi; // Radians

  double m_star; //Global value, but overly important in spot calculations
  double rspot; // Calculated from the equatorial radius and oblateness - unitless
  double Tspot; // In principle can become different - kev

  double redshift_grav; // Equatorial redshift

  bool visible = false;

  std::map<std::string, double> DerivedVals;
  std::map<std::string, std::string> DerivedStrs; // Where possible something should go into the vals.  

  std::map<double, double> BBSpectralFluxMap, HSpectralFluxMap; // Spectral Flux Map
};

struct MainFlags
{
  //Shape function selected and works
  bool shape_model_set = false;
  bool temperature_model_set = true;

  //Basic required inputs set and works
  bool mass_set, req_set, Teq_set, omega_set = false;
  bool inclination_set, distance_set = false;
  bool num_latitude_bins_set, num_longitude_bins_set = false;

  bool model_statistics_header_print = false; 

  bool model_statistics_to_file = false; // Reports to terminal by default
  bool spectra_to_file = false;
  bool grid_to_file = false;
  bool shape_to_file = false;

  bool allInputsSet() {
    return mass_set && req_set && Teq_set && omega_set 
      && inclination_set && distance_set
      && num_latitude_bins_set && num_longitude_bins_set
      && shape_model_set && temperature_model_set;
  }
};

struct ShapeParams
{
  //Shape function related
  int model; 
  std::string model_name = ""; // Name
  std::string shape_library = "shape_functions/lib/PolyOblModelBase"; // File location relative to main running directory
};

//Global star-wise data, needed for the model of the star as a whole.
struct StarParams
{
  double req; // Equatorial Radius. Unitless
  double Teq; // Equatorial Temperature, may be used to set local spot temperatures depending on the temperature model. kev
  double mass; // Overall mass. Unitless
  double omega; // Not modelling differential rotation at the current moment. Unitless

  //Shape function parameters.
  ShapeParams shapeParams;

  //==
  void ReadCMD_Star(int argc, char **argv, MainFlags *MainFlags) {
    // Read in the command line options
    for (int i = 1; i < argc; i++)
    {
      if (argv[i][0] == '-')
      {
        std::string arg = argv[i];
        
        if (arg == "-Sm")
        {
          sscanf(argv[i + 1], "%lf", &mass);
          mass = Units::cgs_to_nounits(mass * Units::MSUN, Units::MASS);
          MainFlags->mass_set = true;
          i++;
        }
        else if (arg == "-Sr")
        {
          sscanf(argv[i + 1], "%lf", &req);
          req = Units::cgs_to_nounits(req * 1.0e5, Units::LENGTH);
          MainFlags->req_set = true;
          i++;
        }
        else if (arg == "-Sf")
        {
          sscanf(argv[i + 1], "%lf", &omega);
          omega = Units::cgs_to_nounits(2.0 * Units::PI * omega, Units::INVTIME);
          MainFlags->omega_set = true;
          i++;
        }
        else if (arg == "-St")
        {
          sscanf(argv[i + 1], "%lf", &Teq);
          MainFlags->Teq_set = true;
          i++;
        }
        else if (arg == "-Ss")
        {
          shapeParams.shape_library = "shape_functions/lib/";
          shapeParams.shape_library += argv[i + 1];
          i++;
        }
        else if (arg == "-Sq")
        {
          sscanf(argv[i + 1], "%d", &shapeParams.model);
          MainFlags->shape_model_set = true;
          i++;
        }
        else if (arg == "-Smn")
        {
          shapeParams.model_name = argv[i + 1];
          MainFlags->shape_model_set = true;
          i++;
        }
      } // end if
    }   // end for
  }
};

//Information on where the observer is
struct ObserverParams
{
  double inclination; // Radians
  double distance; // Unitless

  void ReadCMD_Observer(int argc, char **argv, MainFlags *MainFlags) {
  for (int i = 1; i < argc; i++)
  {
    if (argv[i][0] == '-')
    {
      std::string arg = argv[i];

      if (arg == "-Oi")
      {
        sscanf(argv[i + 1], "%lf", &inclination);
        inclination = fmod(fabs(inclination), 180); //Inclination is only meaningful between 0 and 180 degrees.
        inclination *= Units::PI / 180.0; // To radians.
        MainFlags->inclination_set = true;
        i++;
      }
      else if (arg == "-Od")
      {
        sscanf(argv[i + 1], "%lf", &distance);
        distance = Units::cgs_to_nounits(distance * 100, Units::LENGTH);
        MainFlags->distance_set = true;
        i++;
      }
    } // end if
  }   // end for
}
};

struct GridParams
{
  long num_latitude_bins; // Integer number of bins
  long num_longitude_bins; // Integer number of bins

  void ReadCMD_Grid(int argc, char **argv, MainFlags *MainFlags) {
  for (int i = 1; i < argc; i++)
  {
    if (argv[i][0] == '-')
    {
      std::string arg = argv[i];
      
      if (arg == "-Gn")
      {
        sscanf(argv[i + 1], "%ld", &num_longitude_bins);
        MainFlags->num_longitude_bins_set = true;
        i++;
      }
      else if (arg == "-Gm")
      {
        sscanf(argv[i + 1], "%ld", &num_latitude_bins);
        MainFlags->num_latitude_bins_set = true;
        i++;
      }
    } // end if
  }   // end for
}
};


struct OutputFilesDirectories 
{
  std::string model_statistics_output_file;
  std::string spectra_output_file;
  std::string spot_output_file;
  std::string shape_output_file;

  //==
  void ReadCMD_OutputFiles(int argc, char **argv, MainFlags *MainFlags) {
    // Read in the command line options
    for (int i = 1; i < argc; i++)
    {
      if (argv[i][0] == '-')
      {
        std::string arg = argv[i];
        if (arg == "-Fr")
        {
          model_statistics_output_file = argv[i + 1];
          MainFlags->model_statistics_to_file = true;
          i++;
        }
        else if (arg == "-Fs")
        {
          spectra_output_file = argv[i + 1];
          MainFlags->spectra_to_file = true;
          i++;
        }
        else if (arg == "-Fsh")
        {
          shape_output_file = argv[i + 1];
          MainFlags->shape_to_file = true;
          i++;
        }
        else if (arg == "-Fg")
        {
          spot_output_file = argv[i + 1];
          MainFlags->grid_to_file = true;
          i++;
        }

      }   // end for
    }
  }
};

struct MainParams
{
  //Program boolean flags
  MainFlags mainFlags;

  //Size of the grid
  GridParams gridParams;

  //Stars
  StarParams starParams;

  //Observers
  ObserverParams observerParams;

  //Currently unused
  double quad; // Quadrupole Moment
  double a_kerr;
  double inertia;
  long quadmodel, fdmodel = 0;
  
  OutputFilesDirectories outputFilesDirectories;

  // These two are passed straight to output, useful for labelling
  std::map<std::string, std::string> comments;

  const std::vector<std::string> PRINT_HEADER = {
    "Mass (M_sun)","Equatorial Radius (km)","Rotational Frequency (Hz)","Observer Inclination (Degrees)","Observer Distance (m)","# Latitude Bins","# Longitude Bins","Shape Model Index","Shape Model",
    "Solid Angle (stradians)","Number of Visible Spots","Surface Area (sqkm)",
    "Flux (1keV)",
    "BB Bolometric Flux (Local Area Integral in ergs/cm^2/s)","H Bolometric Flux (Local Area Integral in ergs/cm^2/s)","BB Bolometric Flux (Flux Spectrum Integral in ergs/cm^2/s)","H Bolometric Flux (Flux Spectrum Integral in ergs/cm^2/s)",
    "Steph-Boltzmann BB Flux (ergs/cm^2/s)", "Steph-Boltzmann BB Flux (ergs/cm^2/s) No Doppler",
    "Polar Radius (km)","Polar/Equatorial Radius Ratio",
    "Equatorial Redshift","Polar Redshift",
    "Theoretical Spherical Solid Angle (str)","Approximation Edge-On Ellipsoid Solid Angle (str)",
    "Approximation Edge-On Ellipsoid Solid Angle / Computed Solid Angle",
    "M/R Ratio","Dimensionless Spin Parameter",
    "Furthest Visible Latitude Bin", "Furthest Visible Latitude"
    };

  std::vector<std::string> GetFullHeader() const {
    std::vector<std::string> Run_Output_Header_Final(PRINT_HEADER);
    for (const auto &kv : comments) {
      Run_Output_Header_Final.push_back(kv.first);
    }
    return Run_Output_Header_Final;
  }

  //==

  void ReadCMDArgs(int argc, char **argv)
  {
    // Unpack Parameters
    MainFlags *MainFlags = &mainFlags;

    starParams.ReadCMD_Star(argc, argv, MainFlags);
    observerParams.ReadCMD_Observer(argc, argv, MainFlags);
    gridParams.ReadCMD_Grid(argc, argv, MainFlags);
    outputFilesDirectories.ReadCMD_OutputFiles(argc, argv, MainFlags);

    // Read in the command line options
    for (int i = 1; i < argc; i++)
    {
      if (argv[i][0] == '-')
      {
        std::string arg = argv[i];
        if (arg == "-C")
        {
          comments[argv[i + 1]] = argv[i + 2];
          i+=2;
        }
        if (arg == "-H") {
          std::vector<std::string> Full_Header = GetFullHeader();
          for (std::string item:Full_Header) {
            std::cout << item << ",";
          }
          std::cout << std::endl;
          exit(0);
        }
      } // end if
    }   // end for
  }
};

// Final compilation of all data products.
struct RunOutputs
{
  static inline const std::vector<std::string> ModelStatistics = {"Shape Model", "Omega_s", "Area", "Flux", "BBBoloFluxLocalAreaIntegral",
                                         "HBoloFluxLocalAreaIntegral", "BBBoloFluxSpectrumIntegral",
                                         "HBoloFluxSpectrumIntegral", "BBBoloFluxStephBoltzmann",
                                         "BBBoloFluxStephBoltzmannNoDoppler", "BBBoloFluxStephBoltzmannNoSpin",
                                         "BBBoloFluxStephBoltzmannNoSpinNewtonian",
                                         "ScaledSolidAngle", "UnscaledSolidAngle"};
  
  std::map<double, double> BBSpectralFluxMap, HSpectralFluxMap; // Final spectral flux map for the model.

  std::map<std::string, std::string> outputs; //All model statistics have to be formatted into a string in the end.
  RunOutputs()
  {
    for (std::string value:ModelStatistics) {
      outputs[value] = ""; // Initialize to no string.
    }
  };
};

inline std::string DoubleToString(double value, int precision) {
    std::ostringstream oss;
    oss << std::scientific << std::setprecision(precision) << value;
    return oss.str();
}

#endif // MAINSTRUCTS_H
