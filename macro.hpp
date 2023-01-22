#ifndef MACRO_HPP
#define MACRO_HPP


#include <iostream>
#include <experimental/filesystem>
#include <iomanip>
#include <cmath>
#include <limits>
#include <armadillo>
#include <string>

/**
 * Define boundary conditions to be implemented.
 */
#define PBC
//#define APBC

/**
 * Kibble-Zurek parameter:
 * Choose the parameter to which KZ is applied
*/
//#define KZ_G
//#define KZ_LAMBDA
#define KZ_KAPPA

/**
 * Observables taken into account.
 */
#define OBS_MAG
#define OBS_ADIABATICITY

/**
 * Define useful macros
 */
#define PI (3.1415926535897932384626433)

#endif
