/***************************************************************************************/
/*                                      Hydrogen.h

   Routines for processing Hydrogen

*/
/***************************************************************************************/

int th_index_nsx(double cos_theta, class LightCurve* mexmcc);

void ReadNSXH(class Atmo* atmo);


void NSXHSpectrum(double cos_theta, double logT, double lgrav, class Atmo* atmo, class Spectrum* spectrum);

void CleanNSXH(class Atmo* atmo, class Spectrum* spectrum);
