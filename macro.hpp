#ifndef MACRO_HPP
#define MACRO_HPP


#include <iostream>
#include <experimental/filesystem>
#include <iomanip>
#include <cmath>
#include <limits>
#include <armadillo>
#include <string>
#include <random>

/**
 * Define boundary conditions to be implemented.
 */
#define PBC
//#define APBC

/**
 * If ROUND_TRIP is defined, stop the KZ protocol when the parameters are all the same as the starting ones.
*/
#define ROUND_TRIP
#define N_CYCLES (1000)

/**
 * Kibble-Zurek parameter:
 * Choose the parameter to which KZ is applied
*/
#define KZ_G
//#define KZ_LAMBDA
//#define KZ_KAPPA

/**
 * Observables taken into account.
 */
#define OBS_MAG
#define OBS_ADIABATICITY

/**
 * Define useful macros
 */
#define PI (3.1415926535897932384626433)
#define I (std::complex<double> (0.0, 1.0))

#endif
