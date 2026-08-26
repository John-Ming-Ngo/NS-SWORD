/*
 * Archived, unscaled PVLS implementation.
 *
 * This loadable compatibility library preserves the exact inference entry
 * points used by historical PVLS runs.  All common C++ behaviour remains in
 * Papigkiotis.cpp; only the class symbol, scripts, and reported model name
 * differ.
 */
#define Papigkiotis Papigkiotis_old
#define PAPIGKIOTIS_SURFACE_SCRIPT "shape_functions/dependencies/Papigkiotis_old/ns_radius.py"
#define PAPIGKIOTIS_DERIVATIVE_SCRIPT "shape_functions/dependencies/Papigkiotis_old/ns_log_derivative.py"
#define PAPIGKIOTIS_MODEL_NAME "PVLS_old"
#define PAPIGKIOTIS_ENFORCE_TRAINED_DOMAIN 0
#define PAPIGKIOTIS_ENABLE_CACHE 0
#include "Papigkiotis.cpp"
