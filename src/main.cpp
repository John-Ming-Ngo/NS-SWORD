// Source (C) Coire Cadeau 2007, all rights reserved.
// Modified and expanded 2024 by John Ming Ngo.

#include "MainStructs.h"
#include "MainPrints.h"
#include "../include/utils.h"

#include "../include/interp_functions.h"

// -- Todo: Refactor this at some point
#include <sys/stat.h>   // For mkdir on POSIX
#include <sys/types.h>  // For mode_t
#include <errno.h>
#include <string>
#include <vector>
#include <cstring>      // For strerror
#include <iomanip>
#include <iostream>     // For debug printing

#ifdef _WIN32
  #include <direct.h>   // For _mkdir
  #define MKDIR(name, mode) _mkdir(name)
#else
  #include <unistd.h>   // For access, etc.
  #define MKDIR(name, mode) mkdir(name, mode)
#endif

void InitializeSpectrumMap(std::map<double, double> *spectrumMap, StarParams *StarParams) {
  double T_kev = StarParams->Teq;
  double m_msun = Units::nounits_to_cgs(StarParams->mass, Units::MASS) / Units::MSUN;
  double req_km = Units::nounits_to_cgs(StarParams->req, Units::LENGTH) / 1.0e5;
  InitializeSpectrumMap(spectrumMap, T_kev, m_msun, req_km);
}

void ObtainModel(OblModelBase **model, StarParams *starParams)
{
  ShapeParams *shapeParams = &starParams->shapeParams;
  int modelType = shapeParams->model;
  
  std::string model_name = starParams->shapeParams.model_name;
  bool model_by_name = model_name.length() != 0;

  std::string libName = shapeParams->shape_library;

#if defined(_WIN32) || defined(_WIN64)
  // Windows: Use LoadLibrary and GetProcAddress
  std::string libFile = libName + ".dll";
  HMODULE hLib = LoadLibrary(libFile.c_str());
  if (!hLib)
  {
    throw std::runtime_error("Failed to load DLL: " + libFile);
  }

  typedef OblModelBase *(*CreateModelFunc)(const double &, const double &, const double &, int);
  typedef OblModelBase *(*CreateModelFuncStr)(const double &, const double &, const double &, std::string);

  CreateModelFunc create_model = nullptr;
  CreateModelFuncStr create_model_str = nullptr;

  if (model_by_name) {
    create_model_str = (CreateModelFuncStr)GetProcAddress(hLib, "createOblModelBaseStr");
    if (!create_model_str)
    {
      FreeLibrary(hLib);
      throw std::runtime_error("Failed to load symbol: createOblModelBaseStr");
    }
  }
  else {
    create_model = (CreateModelFunc)GetProcAddress(hLib, "createOblModelBase");
    if (!create_model)
    {
      FreeLibrary(hLib);
      throw std::runtime_error("Failed to load symbol: createOblModelBase");
    }
  }

#elif defined(__APPLE__)
  // macOS: Use dlopen and dlsym, with .dylib extension
  std::string libFile = libName + ".dylib";
  void *handle = dlopen(libFile.c_str(), RTLD_LAZY);
  if (!handle)
  {
    throw std::runtime_error("Failed to load dylib: " + std::string(dlerror()));
  }

  typedef OblModelBase *(*CreateModelFunc)(const double &, const double &, const double &, int);
  typedef OblModelBase *(*CreateModelFuncStr)(const double &, const double &, const double &, std::string);

  CreateModelFunc create_model = nullptr;
  CreateModelFuncStr create_model_str = nullptr;

  if (model_by_name) {
    create_model_str = (CreateModelFuncStr)dlsym(handle, "createOblModelBaseStr");
    const char *dlsym_error = dlerror();
    if (dlsym_error)
    {
      dlclose(handle);
      throw std::runtime_error(
        "Failed to load symbol: createOblModelBaseStr, " + std::string(dlsym_error)
      );
    }
  }
  else {
    create_model = (CreateModelFunc)dlsym(handle, "createOblModelBase");
    const char *dlsym_error = dlerror();
    if (dlsym_error)
    {
      dlclose(handle);
      throw std::runtime_error(
        "Failed to load symbol: createOblModelBase, " + std::string(dlsym_error)
      );
    }
  }

#else
  // Linux: Use dlopen and dlsym, with .so extension
  std::string libFile = libName + ".so";
  void *handle = dlopen(libFile.c_str(), RTLD_LAZY);
  if (!handle)
  {
    throw std::runtime_error("Failed to load shared object: " + std::string(dlerror()));
  }

  typedef OblModelBase *(*CreateModelFunc)(const double &, const double &, const double &, int);
  typedef OblModelBase *(*CreateModelFuncStr)(const double &, const double &, const double &, std::string);

  CreateModelFunc create_model = nullptr;
  CreateModelFuncStr create_model_str = nullptr;

  if (model_by_name) {
    create_model_str = (CreateModelFuncStr)dlsym(handle, "createOblModelBaseStr");
    const char *dlsym_error = dlerror();
    if (dlsym_error)
    {
      dlclose(handle);
      throw std::runtime_error(
        "Failed to load symbol: createOblModelBaseStr, " + std::string(dlsym_error)
      );
    }
  }
  else {
    create_model = (CreateModelFunc)dlsym(handle, "createOblModelBase");
    const char *dlsym_error = dlerror();
    if (dlsym_error)
    {
      dlclose(handle);
      throw std::runtime_error(
        "Failed to load symbol: createOblModelBase, " + std::string(dlsym_error)
      );
    }
  }

#endif

  // Create the model using the correct function pointer
  if (model_by_name) {
    *model = create_model_str(starParams->req, starParams->mass, starParams->omega, model_name);
  }
  else {
    *model = create_model(starParams->req, starParams->mass, starParams->omega, modelType);
  }
}

// Gives the light curve the information it needs other than spot specific info
void SetupLightCurve_Globals(LightCurve *LightCurve, MainParams MainParams)
{
  // Unpack inputs
  MainFlags MainFlags = MainParams.mainFlags;
  StarParams StarParams = MainParams.starParams;
  ObserverParams ObserverParams = MainParams.observerParams;
  GridParams GridParams = MainParams.gridParams;

  LightCurve->para_read_in = false;

  // curve is a structure holding parameters describing the star and emission properties
  LightCurve->para.mass = StarParams.mass;
  LightCurve->para.omega = StarParams.omega;

  LightCurve->numbins = GridParams.num_longitude_bins;
}

void SetupQuadFDModels(MainParams *MainParams)
{
  if (MainParams->quadmodel == 1)
    MainParams->quad = -0.2;
  else
    MainParams->quad = 0.0;

  if (MainParams->fdmodel == 1)
    MainParams->a_kerr = 0.24;
  else
    MainParams->a_kerr = 0.0;

  if (MainParams->fdmodel == 1)
    MainParams->inertia = sqrt(MainParams->starParams.mass / MainParams->starParams.req) * (1.14 - 2.53 * MainParams->starParams.mass / MainParams->starParams.req + 5.6 * pow(MainParams->starParams.mass / MainParams->starParams.req, 2));
  else
    MainParams->inertia = 0;
}

void ComputeBvPsiTable(LightCurve *LightCurve, OblDeflectionTOA *defltoa, double rspot)
{
  double b_mid;
  // Compute b vs psi lookup table good for the specified M/R and mu
  b_mid = LightCurve->defl.b_max * 0.9;
  LightCurve->defl.b_psi[0] = 0.0;
  LightCurve->defl.psi_b[0] = 0.0;
  for (unsigned int i(1); i < NN + 1; i++)
  { // compute table of b vs psi points
    LightCurve->defl.b_psi[i] = b_mid * i / (NN * 1.0);
    LightCurve->defl.psi_b[i] = defltoa->psi_outgoing(LightCurve->defl.b_psi[i],
                                                      rspot, LightCurve->defl.b_max, LightCurve->defl.psi_max, &LightCurve->problem);
  }
  // For arcane reasons, the table is not evenly spaced.
  for (unsigned int i(NN + 1); i < 3 * NN; i++)
  { // compute table of b vs psi points
    LightCurve->defl.b_psi[i] = b_mid + (LightCurve->defl.b_max - b_mid) / 2.0 * (i - NN) / (NN * 1.0);
    LightCurve->defl.psi_b[i] = defltoa->psi_outgoing(LightCurve->defl.b_psi[i],
                                                      rspot, LightCurve->defl.b_max, LightCurve->defl.psi_max, &LightCurve->problem);
  }
  LightCurve->defl.b_psi[3 * NN] = LightCurve->defl.b_max;
  LightCurve->defl.psi_b[3 * NN] = LightCurve->defl.psi_max;
  // Finished computing lookup table
}

double CalculateRedshift(MainParams MainParams, double rspot, double mu)
{
  double P2 = 0; // Currently plays no role beyond disabling something
  double x;
  x = rspot / MainParams.starParams.mass;
  double F1;
  F1 = -5.0 / 8.0 * (x - 1.0) / (x * (x - 2)) * (2.0 + 6.0 * x - 3.0 * x * x) - 15.0 / 16.0 * x * (x - 2) * log(x / (x - 2.0));
  double Sigma;
  Sigma = pow(x, 2) + pow(MainParams.a_kerr * mu, 2);
  double Nsquare;
  Nsquare = (1.0 - 2.0 * x / Sigma) * (1.0 + 2.0 * MainParams.quad * F1 * P2);
  double redshift = 1.0 / sqrt(Nsquare);

  return redshift;
}

void PrintSpectrum(std::map<double, double> spectrum, std::ostream *out, std::string separator = ",")
{
  *out << "Photon Energy (keV), Flux (Ergs/cm^2/s/hz)" << std::endl;
  for (auto &entry : spectrum)
  {
    *out << entry.first << separator << entry.second << std::endl;
  }
}

//==
bool SetSpotCoordinates(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();  // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double dphi = 2 * Units::PI / (nphi * 1.0);    // Step size in phi (azimuthal angle)
  double dtheta = (Units::PI - 1e-6) / (1.0 * ntheta); // Step size in theta // Epsilon of 1e-6 used to avoid zero errors, not a robust method, TODO.

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      (*StarGrid)[i][j].latitude_bins_coord = i;
      (*StarGrid)[i][j].longitude_bins_coord = j;
      
      (*StarGrid)[i][j].theta = i * dtheta + 1e-6 + 0.5 * dtheta; // TODO: If theta > pi/2, do the hemisphere inversion.
      (*StarGrid)[i][j].phi = j * dphi + 1e-6 + 0.5 * dphi;
    }
  }
  return true;
}

bool SetSpotTemperature(std::vector<std::vector<SpotData>> *StarGrid, double Teq) {
  int ntheta = StarGrid->size();  // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      (*StarGrid)[i][j].Tspot = Teq;
    }
  }
  return true;
}

bool SetSpotMass(std::vector<std::vector<SpotData>> *StarGrid, double Mass) {
  int ntheta = StarGrid->size();  // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      (*StarGrid)[i][j].m_star = Mass;
    }
  }
  return true;
}

bool SetSpotRspot(std::vector<std::vector<SpotData>> *StarGrid, OblModelBase *model) {
  int ntheta = StarGrid->size();  // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      double mu_spot = cos((*StarGrid)[i][j].theta);
      spot.rspot = model->R_at_costheta(mu_spot);
      
      spot.DerivedVals["Equatorial Radius"] = model->Req_nounits();
      spot.DerivedVals["Polar Radius"] = model->R_at_costheta(1);
    }
  }
  return true;
}

bool TransferSpotFundimentals(std::vector<std::vector<SpotData>> *StarGrid) {
    int ntheta = StarGrid->size();  // Number of latitude bins
    int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins
    // Transfer spot fundamentals
    for (int i = 0; i < ntheta; i++) {
        for (int j = 0; j < nphi; j++) {
            SpotData &spot = (*StarGrid)[i][j];
            // Transfer values
            spot.DerivedVals["Mass"] = spot.m_star;
            spot.DerivedVals["Radius"] = spot.rspot;
            spot.DerivedVals["Temperature (keV)"] = spot.Tspot;
            spot.DerivedVals["Theta"] = spot.theta;
            spot.DerivedVals["Phi"] = spot.phi;
            spot.DerivedVals["Latitude Bin"] = spot.latitude_bins_coord;
            spot.DerivedVals["Longitude Bin"] = spot.longitude_bins_coord;
            spot.DerivedVals["Visible?"] = spot.visible; // This doesn't work when it's executed before the function that puts this in! Hmm.
            spot.DerivedVals["Redshift Factor"] = spot.redshift_grav; // This doesn't work when redshift is being calculated after... Well that's annoying. Time to fix.
        }
    }
    return true;
}

//==

bool CalculateSpotCosG(std::vector<std::vector<SpotData>> *StarGrid, OblModelBase *model) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      double mu_spot = cos((*StarGrid)[i][j].theta);
      spot.DerivedVals["gamma"] =  acos(model->cos_gamma(mu_spot));
      spot.DerivedVals["cosg"] =  model->cos_gamma(mu_spot);
      spot.DerivedVals["sing"] =  sin(spot.DerivedVals["gamma"]);
    }
  }
  return true;
}

bool CalculateSpot_dS(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();  // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double dphi = 2 * Units::PI / (nphi * 1.0);    // Step size in phi (azimuthal angle)
  double dtheta = (Units::PI - 1e-6) / (1.0 * ntheta); // Step size in theta // Epsilon of 1e-6 used to avoid zero errors, not a robust method, TODO.

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      (*StarGrid)[i][j].DerivedVals["dS"] = (sin((*StarGrid)[i][j].theta) * (dtheta * dphi)) / (*StarGrid)[i][j].DerivedVals["cosg"];
    }
  }
  return true;  
}

bool CalculateSpotLightCurve(std::vector<std::vector<SpotData>> *StarGrid, LightCurve *LightCurve, OblDeflectionTOA *defltoa, double inclination) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double dphi = 2 * Units::PI / (nphi * 1.0);    // Step size in phi (azimuthal angle)
  double dtheta = (Units::PI - 1e-6) / (1.0 * ntheta); // Step size in theta //Epsilon of 1e-6 used to avoid zero errors, not a robust method, TODO.

  for (int i = 0; i < ntheta; i++) {
    double theta = i * dtheta + 1e-6 + 0.5 * dtheta;
    double mu = cos(theta);

    LightCurve->para.theta = theta;
    LightCurve->para.dS = (*StarGrid)[i][0].DerivedVals["dS"];
    LightCurve->para.cosgamma = (*StarGrid)[i][0].DerivedVals["cosg"];
    LightCurve->para.radius = (*StarGrid)[i][0].rspot;

    //std::cout << inclination << std::endl;
    LightCurve->para.incl = inclination;

    // Calculate deflection values
    LightCurve->defl.b_max = defltoa->bmax_outgoing((*StarGrid)[i][0].rspot);
    LightCurve->defl.psi_max = defltoa->psi_max_outgoing(LightCurve->defl.b_max, (*StarGrid)[i][0].rspot, &LightCurve->problem);

    // Compute table of b vs psi points
    ComputeBvPsiTable(LightCurve, defltoa, (*StarGrid)[i][0].rspot); //Read only on defltoa

    // Compute angles for this latitude point
    *LightCurve = ComputeAngles(LightCurve, defltoa); //Read only on defltoa
    for (int j = 0; j < nphi; j++) {
      SpotData &Spot = (*StarGrid)[i][j];
      Spot.DerivedVals["dOmega_s"] = LightCurve->dOmega_s[j];
      Spot.DerivedVals["eta"] = LightCurve->eta[j];
      
      Spot.DerivedVals["beta"] = LightCurve->beta[j];
      Spot.DerivedVals["sinbeta"] = LightCurve->sinbeta[j];
      Spot.DerivedVals["cosbeta"] = LightCurve->cosbeta[j];

      Spot.DerivedVals["alpha"] = LightCurve->alpha[j];
      Spot.DerivedVals["sinalpha"] = LightCurve->sinalpha[j];
      Spot.DerivedVals["cosalpha"] = LightCurve->cosalpha[j];

      Spot.visible = LightCurve->visible[j];
      Spot.DerivedVals["Visible?"] = Spot.visible; // Todo: This is a dissatisfying patch job.

      Spot.DerivedVals["t"] = LightCurve->t[j];
      Spot.DerivedVals["t_o"] = LightCurve->t_o[j];
      Spot.DerivedVals["cosxi"] = LightCurve->cosxi[j];

      Spot.DerivedVals["Impact Parameter b"] = LightCurve->b[j];
      Spot.DerivedVals["psi"] = LightCurve->psi[j];
      Spot.DerivedVals["cosdelta"] = LightCurve->cosdelta[j];
      Spot.DerivedVals["dcosalpha_dcospsi"] = LightCurve->dcosalpha_dcospsi[j];
    }
  }
  return true;
}

bool CalculateSpotLightCurve(std::vector<std::vector<SpotData>> *StarGrid, LightCurve *LightCurve, OblDeflectionTOA *defltoa, double inclination, double req) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double dphi = 2 * Units::PI / (nphi * 1.0);    // Step size in phi (azimuthal angle)
  double dtheta = (Units::PI - 1e-6) / (1.0 * ntheta); // Step size in theta //Epsilon of 1e-6 used to avoid zero errors, not a robust method, TODO.

  for (int i = 0; i < ntheta; i++) {
    double theta = i * dtheta + 1e-6 + 0.5 * dtheta;
    double mu = cos(theta);

    LightCurve->para.theta = theta;
    LightCurve->para.dS = (*StarGrid)[i][0].DerivedVals["dS"];
    LightCurve->para.cosgamma = (*StarGrid)[i][0].DerivedVals["cosg"];
    LightCurve->para.radius = req;

    //std::cout << inclination << std::endl;
    LightCurve->para.incl = inclination;

    // Calculate deflection values
    LightCurve->defl.b_max = defltoa->bmax_outgoing((*StarGrid)[i][0].rspot);
    LightCurve->defl.psi_max = defltoa->psi_max_outgoing(LightCurve->defl.b_max, (*StarGrid)[i][0].rspot, &LightCurve->problem);

    // Compute table of b vs psi points
    ComputeBvPsiTable(LightCurve, defltoa, (*StarGrid)[i][0].rspot); //Read only on defltoa

    // Compute angles for this latitude point
    *LightCurve = ComputeAngles(LightCurve, defltoa); //Read only on defltoa
    for (int j = 0; j < nphi; j++) {
      SpotData &Spot = (*StarGrid)[i][j];
      Spot.DerivedVals["dOmega_s"] = LightCurve->dOmega_s[j];
      Spot.DerivedVals["eta"] = LightCurve->eta[j];
      Spot.DerivedVals["cosbeta"] = LightCurve->cosbeta[j];

      Spot.visible = LightCurve->visible[j];
      Spot.DerivedVals["Visible?"] = Spot.visible; // Todo: This is a dissatisfying patch job.

      Spot.DerivedVals["t"] = LightCurve->t[j];
      Spot.DerivedVals["t_o"] = LightCurve->t_o[j];
      Spot.DerivedVals["cosxi"] = LightCurve->cosxi[j];

      Spot.DerivedVals["Impact Parameter b"] = LightCurve->b[j];
      Spot.DerivedVals["psi"] = LightCurve->psi[j];
      Spot.DerivedVals["cosdelta"] = LightCurve->cosdelta[j];
      Spot.DerivedVals["dcosalpha_dcospsi"] = LightCurve->dcosalpha_dcospsi[j];
    }
  }
  return true;
}

bool CalculateSpotRedshift(std::vector<std::vector<SpotData>> *StarGrid, MainParams MainParams) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      (*StarGrid)[i][j].redshift_grav = CalculateRedshift(MainParams, (*StarGrid)[i][j].rspot, cos((*StarGrid)[i][j].theta));
    }
  }
  return true;
}

bool CalculateSpotArea(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      (*StarGrid)[i][j].DerivedVals["dArea"] = (*StarGrid)[i][j].DerivedVals["dS"] * pow(Units::nounits_to_cgs((*StarGrid)[i][j].rspot * 1.0e-5, Units::LENGTH), 2);
      /*
      double theta = (*StarGrid)[i][j].theta;
      double sinTheta = sin(theta);
      double cosTheta = cos(theta);

      double equatorial_radius = Units::nounits_to_cgs((*StarGrid)[i][j].DerivedVals["Equatorial Radius"] * 1.0e-5, Units::LENGTH);
      double polar_radius = Units::nounits_to_cgs((*StarGrid)[i][j].DerivedVals["Polar Radius"] * 1.0e-5, Units::LENGTH);

      (*StarGrid)[i][j].DerivedVals["dArea"] = (*StarGrid)[i][j].DerivedVals["dS"] 
      * equatorial_radius
      * sqrt(pow(equatorial_radius * cosTheta, 2) + pow(polar_radius * sinTheta, 2));
      */
    }
  }
  return true;
}

bool CalculateSpotScaledSolidAngle(std::vector<std::vector<SpotData>> *StarGrid, double distance) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &Spot = (*StarGrid)[i][j];
      if (Spot.DerivedVals["dOmega_s"] != 0.0) {
        Spot.DerivedVals["dOmega"] = Spot.DerivedVals["dOmega_s"] * pow(Spot.rspot / distance, 2);
        /*
        double theta = (*StarGrid)[i][j].theta;
        double sinTheta = sin(theta);
        double cosTheta = cos(theta);

        double equatorial_radius = (*StarGrid)[i][j].DerivedVals["Equatorial Radius"];
        double polar_radius = (*StarGrid)[i][j].DerivedVals["Polar Radius"];

        double area_scaling = equatorial_radius * sqrt(pow(equatorial_radius * cosTheta, 2) + pow(polar_radius * sinTheta, 2));
        Spot.DerivedVals["dOmega"] = Spot.DerivedVals["dOmega_s"] * pow(1/ distance, 2) * area_scaling;
        */
      }
      else {
        (*StarGrid)[i][j].DerivedVals["dOmega"] = 0.0;
      }
    }
  }
  return true;
}

//==

bool CalculateSpotFlux_Legacy(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      if ((*StarGrid)[i][j].DerivedVals["dOmega_s"] != 0.0) {
        (*StarGrid)[i][j].DerivedVals["dBoloFlux_mu"] =  (*StarGrid)[i][j].DerivedVals["dOmega"] * pow((*StarGrid)[i][j].DerivedVals["eta"], 4) / Units::PI * pow((*StarGrid)[i][j].redshift_grav, -4);
        (*StarGrid)[i][j].DerivedVals["dFlux_mu"] = (*StarGrid)[i][j].DerivedVals["dOmega"] * pow((*StarGrid)[i][j].DerivedVals["eta"], 3) * BlackBody((*StarGrid)[i][j].Tspot, 1.0 * (*StarGrid)[i][j].redshift_grav / (*StarGrid)[i][j].DerivedVals["eta"]) * pow(1 / (*StarGrid)[i][j].redshift_grav, 3) * (1.0 / (1.0 * Units::H_PLANCK));
      }
      else {
        (*StarGrid)[i][j].DerivedVals["dBoloFlux_mu"] = 0.0;
        (*StarGrid)[i][j].DerivedVals["dFlux_mu"] = 0.0;
      }
    }
  }
  return true;
}

//==

bool InitializeStarGridSpectra(std::vector<std::vector<SpotData>> *StarGrid, StarParams *StarParams) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      InitializeSpectrumMap(&((*StarGrid)[i][j].BBSpectralFluxMap), StarParams);
      InitializeSpectrumMap(&((*StarGrid)[i][j].HSpectralFluxMap), StarParams);
    }
  }
  return true;
}

bool CalculateSpotBBSpectra(std::vector<std::vector<SpotData>> *StarGrid, double mass) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      if ((*StarGrid)[i][j].DerivedVals["dOmega_s"] != 0.0)
      {
        std::map<double, double> BBSpectrum;
        BBSpectrum.clear();
        GetBBSpectrumMap(&BBSpectrum,
                         (*StarGrid)[i][j].DerivedVals["cosbeta"],
                         (*StarGrid)[i][j].Tspot,
                         Units::nounits_to_cgs(mass, Units::MASS) / Units::MSUN,
                         Units::nounits_to_cgs((*StarGrid)[i][j].rspot, Units::LENGTH) / 1.0e5
        );
        spectrum_redshift(&BBSpectrum, (*StarGrid)[i][j].DerivedVals["eta"] / (*StarGrid)[i][j].redshift_grav);

        (*StarGrid)[i][j].DerivedVals["Local Blackbody Bolometric Flux Spectrum Integral Value"] =
            spectrum_integrate(BBSpectrum, 1) * (*StarGrid)[i][j].DerivedVals["dOmega"] * Units::EV_TO_HZ * 1000;             // Maybe redshift light curve expansion issue related?
                                                     // Are angles somehow being treated differently between the two? how?

        spectrum_multiply(&BBSpectrum, (*StarGrid)[i][j].DerivedVals["dOmega"] * Units::EV_TO_HZ * 1000);
        spectrum_interpolate_add(&(*StarGrid)[i][j].BBSpectralFluxMap, &BBSpectrum, 3);
      }
    }
  }
  return true;
}

bool CalculateSpotHSpectra(std::vector<std::vector<SpotData>> *StarGrid, double mass) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      if (spot.DerivedVals["dOmega_s"] != 0.0)
      {
        std::map<double, double> HSpectrum;
        GetHSpectrumMap(&HSpectrum,
                        spot.DerivedVals["cosbeta"],
                        spot.Tspot,
                        Units::nounits_to_cgs(mass, Units::MASS) / Units::MSUN,
                        Units::nounits_to_cgs(spot.rspot, Units::LENGTH) / 1.0e5);

        spectrum_redshift(&HSpectrum, spot.DerivedVals["eta"] / spot.redshift_grav);

        spot.DerivedVals["Local Hydrogen Bolometric Flux Spectrum Integral Value"] =
            spectrum_integrate(HSpectrum, 1) * spot.DerivedVals["dOmega"] * Units::EV_TO_HZ * 1000;

        spectrum_multiply(&HSpectrum, spot.DerivedVals["dOmega"] * Units::EV_TO_HZ * 1000);
        spectrum_interpolate_add(&spot.HSpectralFluxMap, &HSpectrum, 3);
      }
    }
  }
  return true;
}

//==
bool CalculateSpotSBFlux(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      spot.DerivedVals["Steph-Boltzmann BB Flux (ergs/cm^2/s)"] = StephBoltzmann(spot.Tspot * Units::EV_TO_K * 1000, 1.0) 
        * pow(spot.DerivedVals["eta"], 4) / pow(spot.redshift_grav, 4) 
        / Units::PI 
        * spot.DerivedVals["dOmega"]; //Has two powers of redshift in there
    }
  }
  return true;
}

bool CalculateSpotSBFluxNoDoppler(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      spot.DerivedVals["Steph-Boltzmann BB Flux (ergs/cm^2/s) No Doppler"] = StephBoltzmann(spot.Tspot * Units::EV_TO_K * 1000, 1.0) 
        / pow(spot.redshift_grav, 4) 
        / Units::PI 
        * spot.DerivedVals["dOmega"];
    }
  }
  return true;
}

//==
double SumSolidAngles(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double SolidAngle = 0.0;
  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      SolidAngle += spot.DerivedVals["dOmega"];
    }
  }
  return SolidAngle;
}

double SumSurfaceAreas(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double SurfaceArea = 0.0;
  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      SurfaceArea += spot.DerivedVals["dArea"];
    }
  }
  return SurfaceArea;
}

double SumFluxLegacy(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double Flux_1keV = 0.0;
  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      Flux_1keV += spot.DerivedVals["dFlux_mu"];
    }
  }
  return Flux_1keV;
}

double SumLocalBlackBodyFlux(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double LocalBlackbodyFlux = 0.0;
  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      LocalBlackbodyFlux += spot.DerivedVals["Local Blackbody Bolometric Flux Spectrum Integral Value"];
    }
  }
  return LocalBlackbodyFlux;
}

double SumLocalHydrogenFlux(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double LocalBlackbodyFlux = 0.0;
  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      LocalBlackbodyFlux += spot.DerivedVals["Local Hydrogen Bolometric Flux Spectrum Integral Value"];
    }
  }
  return LocalBlackbodyFlux;
}

void SumBBSpectra(std::vector<std::vector<SpotData>> *StarGrid, std::map<double, double> *BBSpectralFluxMap) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  std::map<double, double> &BBSpectra = *BBSpectralFluxMap;

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      spectrum_interpolate_add(&BBSpectra, &spot.BBSpectralFluxMap, 3);
    }
  }
}

void SumHSpectra(std::vector<std::vector<SpotData>> *StarGrid, std::map<double, double> *HSpectralFluxMap) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  std::map<double, double> &Spectra = *HSpectralFluxMap;

  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      spectrum_interpolate_add(&Spectra, &spot.HSpectralFluxMap, 3);
    }
  }
}

double SumStephBoltz(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double StephBoltzFlux = 0.0;
  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      StephBoltzFlux += spot.DerivedVals["Steph-Boltzmann BB Flux (ergs/cm^2/s)"];
    }
  }
  return StephBoltzFlux;
}

double SumStephBoltzNoDoppler(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double StephBoltzFluxNoDoppler = 0.0;
  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      StephBoltzFluxNoDoppler += spot.DerivedVals["Steph-Boltzmann BB Flux (ergs/cm^2/s) No Doppler"];
    }
  }
  return StephBoltzFluxNoDoppler;
}

double SumNumVisibleSpots(std::vector<std::vector<SpotData>> *StarGrid) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  double NumVisible = 0.0;
  for (int i = 0; i < ntheta; i++) {
    for (int j = 0; j < nphi; j++) {
      SpotData &spot = (*StarGrid)[i][j];
      if (spot.visible) NumVisible += 1.0;
    }
  }
  return NumVisible;
}

double angular_distance(double theta1, double phi1, double theta2, double phi2) {
    return acos(sin(theta1) * sin(theta2) * cos(phi1 - phi2) + cos(theta1) * cos(theta2));
}

double FindFurthestVisible(std::vector<std::vector<SpotData>> *StarGrid, double inclination) {
  int ntheta = StarGrid->size();      // Number of latitude bins
  int nphi = ntheta > 0 ? (*StarGrid)[0].size() : 0; // Number of longitude bins

  int behind_the_star = nphi/2;
  int behind_the_star_even_bins_fix = (nphi%2 == 0) ? nphi/2 + 1 : nphi/2; // Might not be necessary if all we care about is the theta angle.

  double max_angle = 0;
  double coord = 0;
  for (int i = 0; i < ntheta; i++) {
    double angle = angular_distance((*StarGrid)[i][behind_the_star].theta, (*StarGrid)[i][behind_the_star].phi, inclination, 0);
    if (((*StarGrid)[i][behind_the_star].visible > 0) && (angle > max_angle)) {
      max_angle = angle;
      coord = i;
    }
  }
  return coord;  
}
//==

template <typename T>
std::vector<std::string> getKeysLoop(const std::map<std::string, T>& myMap)
{
    std::vector<std::string> keys;
    keys.reserve(myMap.size());
    for (const auto& kv : myMap)
    {
        keys.push_back(kv.first);
    }
    return keys;
}

//==

void write_shape_file(const MainParams& MainParams, OblModelBase* model) {
    std::string outFile = MainParams.outputFilesDirectories.shape_output_file;
    std::string dirPath = outFile.substr(0, outFile.find_last_of("/\\"));

    if (!create_directories(dirPath)) {
        std::cerr << "Failed to create directories for " << dirPath << std::endl;
    }

    std::ofstream out;
    out.open(outFile, std::ios::trunc);
    out << std::setprecision(17);

    int nmu = MainParams.gridParams.num_latitude_bins;
    double dtheta = (Units::PI - 1e-6) / (1.0 * nmu);
    out << "Mu, Rspot, Rspot_cgs (km), Req, Req_cgs (km), Ratio, dRdT, dRdT (km)\n";
    double previous_mu = 1.0;
    for (int j = 0; j < nmu; j++) {
        double theta = j * dtheta + 1e-6 + 0.5 * dtheta;
        double mu = cos(theta);

        // A midpoint grid does not include the equator. Record the model's
        // exact boundary value when the traversal crosses from north to south.
        if (previous_mu > 0.0 && mu < 0.0) {
            double equator_rspot = model->R_at_costheta(0.0);
            double equator_rspot_cgs = Units::nounits_to_cgs(equator_rspot, Units::LENGTH) / 1.0e5;
            double equator_req = model->Req_nounits();
            double equator_req_cgs = Units::nounits_to_cgs(equator_req, Units::LENGTH) / 1.0e5;
            double equator_drdt = model->Dtheta_R(0.0);
            double equator_drdt_cgs = Units::nounits_to_cgs(equator_drdt, Units::LENGTH) / 1.0e5;
            out << 0.0 << "," << equator_rspot << "," << equator_rspot_cgs << "," << equator_req << "," << equator_req_cgs
                << "," << equator_rspot/equator_req << "," << equator_drdt << "," << equator_drdt_cgs << "\n";
        }

        double rspot = model->R_at_costheta(mu);
        double rspot_cgs = Units::nounits_to_cgs(rspot, Units::LENGTH) / 1.0e5;
        double req = model->Req_nounits();
        double req_cgs = Units::nounits_to_cgs(req, Units::LENGTH) / 1.0e5;
        double drdt = model->Dtheta_R(mu);
        double drdt_cgs = Units::nounits_to_cgs(drdt, Units::LENGTH) / 1.0e5;
        out << mu << "," << rspot << "," << rspot_cgs << "," << req << "," << req_cgs
            << "," << rspot/model->Req_nounits() << "," << drdt << "," << drdt_cgs << "\n";
        previous_mu = mu;
    }
    out.close();
}

void write_shape_file_only(const MainParams& MainParams) {
    OblModelBase *model = nullptr;
    ObtainModel(&model, const_cast<StarParams*>(&MainParams.starParams)); // Non-const API
    write_shape_file(MainParams, model);
    delete model;
}

void write_model_statistics(MainParams& MainParams, RunOutputs& RunOutputs) {
    std::string outFile = MainParams.outputFilesDirectories.model_statistics_output_file;
    std::string dirPath = outFile.substr(0, outFile.find_last_of("/\\"));

    if (!create_directories(dirPath)) {
        std::cerr << "Failed to create directories for " << dirPath << std::endl;
    }
    std::ofstream out;
    out.open(outFile, std::ios::trunc);
    PrintOutput(RunOutputs.outputs, MainParams.GetFullHeader(), ",", out);
    out.close();
}

void write_spectra(const MainParams& MainParams, const RunOutputs& RunOutputs) {
    std::string outFile = MainParams.outputFilesDirectories.spectra_output_file;
    std::string dirPath = outFile.substr(0, outFile.find_last_of("/\\"));
    if (!create_directories(dirPath)) {
        std::cerr << "Failed to create directories for " << dirPath << std::endl;
    }

    std::ofstream bbout;
    bbout.open(outFile + "BB.csv", std::ios::trunc);
    PrintSpectrum(RunOutputs.BBSpectralFluxMap, &bbout);
    bbout.close();

    std::ofstream hout;
    hout.open(outFile + "H.csv", std::ios::trunc);
    PrintSpectrum(RunOutputs.HSpectralFluxMap, &hout);
    hout.close();
}

void write_spot_grid(const MainParams& MainParams, std::vector<std::vector<SpotData>>& StarGrid) {
    std::string outFile = MainParams.outputFilesDirectories.spot_output_file;
    std::string dirPath = outFile.substr(0, outFile.find_last_of("/\\"));
    if (!create_directories(dirPath)) {
        std::cerr << "Failed to create directories for " << dirPath << std::endl;
    }

    std::ofstream out;
    out.open(outFile, std::ios::trunc);
    out << std::setprecision(9);
    std::vector<std::string> keys = getKeysLoop(StarGrid[0][0].DerivedVals);
    PrintSpotData(&StarGrid, keys, ",", out);
    out.close();
}



int main(int argc, char **argv)
try
{
  // == Initial Setup == 
  MainParams MainParams;   // Program inputs
  
  // Read command-line arguments
  MainParams.ReadCMDArgs(argc, argv);
  
  // Ensure inputs are correct
  if (!MainParams.mainFlags.allInputsSet()) PrintMainFlags(&MainParams.mainFlags);

  // Early exit: only shape output wanted
  if (
      MainParams.mainFlags.shape_to_file &&
      !MainParams.mainFlags.model_statistics_to_file &&
      !MainParams.mainFlags.spectra_to_file &&
      !MainParams.mainFlags.grid_to_file
    ) {
        write_shape_file_only(MainParams);
        return 0;
  }

  // Initialize Global Data Containers (ex. what bins they should be recording for)
  SetupQuadFDModels(&MainParams); // Todo: This is a terrible function, and should be replaced. Also doesn't do anything relevant yet.

  OblModelBase *model = nullptr;
  ObtainModel(&model, &(MainParams.starParams));
  OblDeflectionTOA *defltoa = new OblDeflectionTOA(model, MainParams.starParams.mass);

  //Latitude specific data containers
  LightCurve LightCurve;
  SetupLightCurve_Globals(&LightCurve, MainParams); 
  // Gives it star and observer info; does not give spot-specific data as that'll change.

  // == Core Calculations == 
  std::vector<std::vector<SpotData>> StarGrid(MainParams.gridParams.num_latitude_bins, std::vector<SpotData>(MainParams.gridParams.num_longitude_bins));

  //==
  SetSpotCoordinates(&StarGrid);
  SetSpotMass(&StarGrid, MainParams.starParams.mass);
  SetSpotTemperature(&StarGrid, MainParams.starParams.Teq);
  
  //==
  SetSpotRspot(&StarGrid, model);

  CalculateSpotCosG(&StarGrid, model);
  CalculateSpot_dS(&StarGrid);
  CalculateSpotRedshift(&StarGrid, MainParams);

  //==
  TransferSpotFundimentals(&StarGrid);

  //==
  CalculateSpotLightCurve(&StarGrid, &LightCurve, defltoa, MainParams.observerParams.inclination);
  //CalculateSpotLightCurve(&StarGrid, &LightCurve, defltoa, MainParams.observerParams.inclination, MainParams.starParams.req);

  //==
  CalculateSpotArea(&StarGrid);
  CalculateSpotScaledSolidAngle(&StarGrid, MainParams.observerParams.distance);
  CalculateSpotFlux_Legacy(&StarGrid);
  
  //==
  InitializeStarGridSpectra(&StarGrid, &MainParams.starParams);
  CalculateSpotBBSpectra(&StarGrid, MainParams.starParams.mass);
  CalculateSpotHSpectra(&StarGrid, MainParams.starParams.mass);

  //==
  CalculateSpotSBFlux(&StarGrid);
  CalculateSpotSBFluxNoDoppler(&StarGrid);
  
  // Output model statistics and/or spectra maps.
  RunOutputs RunOutputs; // Program outputs

  //==
  InitializeSpectrumMap(&RunOutputs.BBSpectralFluxMap, &MainParams.starParams);
  InitializeSpectrumMap(&RunOutputs.HSpectralFluxMap, &MainParams.starParams);
  SumBBSpectra(&StarGrid, &RunOutputs.BBSpectralFluxMap);
  SumHSpectra(&StarGrid, &RunOutputs.HSpectralFluxMap);

  /**
   * What gets printed out is handled by the Print Header in MainStructs.h!
   */
  //==
  int OUTPUT_PRECISION = 17;
  StarParams StarParams = MainParams.starParams;
  ShapeParams ShapeParams = StarParams.shapeParams;
  ObserverParams ObserverParams = MainParams.observerParams;
  GridParams GridParams = MainParams.gridParams;
  
  RunOutputs.outputs["Mass (M_sun)"] = DoubleToString(Units::nounits_to_cgs(StarParams.mass, Units::MASS) / Units::MSUN, OUTPUT_PRECISION);
  RunOutputs.outputs["Equatorial Radius (km)"] = DoubleToString(Units::nounits_to_cgs(StarParams.req, Units::LENGTH) / 1.0e5, OUTPUT_PRECISION);
  RunOutputs.outputs["Rotational Frequency (Hz)"] = DoubleToString(Units::nounits_to_cgs(StarParams.omega, Units::INVTIME) / (2.0 * Units::PI), OUTPUT_PRECISION);
  RunOutputs.outputs["Observer Inclination (Degrees)"] = DoubleToString(ObserverParams.inclination/(Units::PI / 180.0), OUTPUT_PRECISION);
  RunOutputs.outputs["Observer Distance (m)"] = DoubleToString(Units::nounits_to_cgs(ObserverParams.distance, Units::LENGTH) / 100, OUTPUT_PRECISION);
  RunOutputs.outputs["# Latitude Bins"] = std::to_string(GridParams.num_latitude_bins);
  RunOutputs.outputs["# Longitude Bins"] = std::to_string(GridParams.num_longitude_bins);
  RunOutputs.outputs["Shape Model Index"] = std::to_string(ShapeParams.model);
  model->get_model_name(&RunOutputs.outputs["Shape Model"]);

  //==
  RunOutputs.outputs["Solid Angle (stradians)"] = DoubleToString(SumSolidAngles(&StarGrid), OUTPUT_PRECISION);
  RunOutputs.outputs["Number of Visible Spots"] = DoubleToString(SumNumVisibleSpots(&StarGrid), OUTPUT_PRECISION);
  RunOutputs.outputs["Surface Area (sqkm)"] = DoubleToString(SumSurfaceAreas(&StarGrid), OUTPUT_PRECISION);

  //==
  RunOutputs.outputs["Flux (1keV)"] = DoubleToString(SumFluxLegacy(&StarGrid), OUTPUT_PRECISION);

  //==
  RunOutputs.outputs["BB Bolometric Flux (Local Area Integral in ergs/cm^2/s)"] = DoubleToString(SumLocalBlackBodyFlux(&StarGrid), OUTPUT_PRECISION);
  RunOutputs.outputs["H Bolometric Flux (Local Area Integral in ergs/cm^2/s)"] = DoubleToString(SumLocalHydrogenFlux(&StarGrid), OUTPUT_PRECISION);
  RunOutputs.outputs["BB Bolometric Flux (Flux Spectrum Integral in ergs/cm^2/s)"] = DoubleToString(spectrum_integrate(RunOutputs.BBSpectralFluxMap), OUTPUT_PRECISION);
  RunOutputs.outputs["H Bolometric Flux (Flux Spectrum Integral in ergs/cm^2/s)"] = DoubleToString(spectrum_integrate(RunOutputs.HSpectralFluxMap), OUTPUT_PRECISION);

  RunOutputs.outputs["Steph-Boltzmann BB Flux (ergs/cm^2/s)"] = DoubleToString(SumStephBoltz(&StarGrid), OUTPUT_PRECISION);
  RunOutputs.outputs["Steph-Boltzmann BB Flux (ergs/cm^2/s) No Doppler"] = DoubleToString(SumStephBoltzNoDoppler(&StarGrid), OUTPUT_PRECISION);

  //==
  RunOutputs.outputs["Polar Radius (km)"] = DoubleToString(Units::nounits_to_cgs(model->R_at_costheta(cos(0)), Units::LENGTH) / 1.0e5, OUTPUT_PRECISION);
  RunOutputs.outputs["Polar/Equatorial Radius Ratio"] = DoubleToString(model->R_at_costheta(1)/model->R_at_costheta(0), OUTPUT_PRECISION);

  //==
  double equatorial_radius = model->R_at_costheta(0);
  double polar_radius = model->R_at_costheta(1);
  double equatorial_redshift = CalculateRedshift(MainParams, model->R_at_costheta(0), 0);
  double polar_redshift = CalculateRedshift(MainParams, model->R_at_costheta(1), 1);

  double spherical_solid_angle = Units::PI * pow(equatorial_radius, 2) * pow(equatorial_redshift, 2);
  double elipsoid_solid_angle = Units::PI * equatorial_radius * polar_radius * (equatorial_redshift * polar_redshift);

  double distance_scaling = pow(ObserverParams.distance, 2);

  RunOutputs.outputs["Equatorial Redshift"] = DoubleToString(equatorial_redshift, OUTPUT_PRECISION);
  RunOutputs.outputs["Polar Redshift"] = DoubleToString(polar_redshift, OUTPUT_PRECISION);
  
  RunOutputs.outputs["Theoretical Spherical Solid Angle (str)"] = DoubleToString(spherical_solid_angle/distance_scaling, OUTPUT_PRECISION);
  RunOutputs.outputs["Approximation Edge-On Ellipsoid Solid Angle (str)"] = DoubleToString(elipsoid_solid_angle/distance_scaling, OUTPUT_PRECISION);
  RunOutputs.outputs["Approximation Edge-On Ellipsoid Solid Angle / Computed Solid Angle"] = DoubleToString((elipsoid_solid_angle/distance_scaling)/SumSolidAngles(&StarGrid), OUTPUT_PRECISION);

  //==
  RunOutputs.outputs["M/R Ratio"] = DoubleToString(StarParams.mass/StarParams.req, OUTPUT_PRECISION);
  RunOutputs.outputs["Dimensionless Spin Parameter"] = DoubleToString(double(pow(StarParams.omega,2.0) * pow(StarParams.req,3.0) /StarParams.mass), OUTPUT_PRECISION);

  //==
  RunOutputs.outputs["Furthest Visible Latitude Bin"] = DoubleToString(FindFurthestVisible(&StarGrid, ObserverParams.inclination), OUTPUT_PRECISION);
  RunOutputs.outputs["Furthest Visible Latitude"] = DoubleToString(FindFurthestVisible(&StarGrid, ObserverParams.inclination) * 180.0 / GridParams.num_latitude_bins, OUTPUT_PRECISION);

  //== Todo: This sucks, should do so in a more intuitive manner.
  for (const auto &kv : MainParams.comments) {
    RunOutputs.outputs[kv.first] = kv.second;
  }

  //==
  if (MainParams.mainFlags.model_statistics_to_file) write_model_statistics(MainParams, RunOutputs);
  else PrintOutput(RunOutputs.outputs, MainParams.GetFullHeader());

  if (MainParams.mainFlags.spectra_to_file) write_spectra(MainParams, RunOutputs);

  //std::cout << "Is supposed to output data files: " << MainParams.mainFlags.grid_to_file << std::endl;
  if (MainParams.mainFlags.grid_to_file) write_spot_grid(MainParams,StarGrid);
  if (MainParams.mainFlags.shape_to_file) write_shape_file(MainParams, model);
  delete defltoa;
  delete model;
  return 0;
}
catch (std::exception &e)
{
  std::cerr << "Top-level exception caught in application: " << e.what() << std::endl;
  return -1;
}
