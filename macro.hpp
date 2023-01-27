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
 * From here we can choose between different schemes for time evolution.
 * -  KZ_PROTOCOL : KZ of one parameter of the hamiltonian
 * - MEASUREMENT_PROTOCOL : The dynamic is triggered by monitoring continuously the central spin system.
*/
//#define KZ_PROTOCOL
#define ROUND_TRIP
#define N_CYCLES (1000)

#define MEASUREMENT_PROTOCOL

/**
 * Kibble-Zurek parameter:
 * Choose the parameter to which KZ is applied
*/
#define KZ_G
//#define KZ_LAMBDA
//#define KZ_KAPPA

/**
 * Measurement dyanmics protocol used
*/
//#define MEASURE_SIGMAX
#define MEASURE_SIGMAZ
//#define MEASURE_RANDOM


/**
 * Observables taken into account.
 */
#define OBS_MAG
//#define OBS_ADIABATICITY

/**
 * Define useful macros
 */
#define PI (3.1415926535897932384626433)
#define I (std::complex<double> (0.0, 1.0))
#define EPSILON (1e-15)

#endif
