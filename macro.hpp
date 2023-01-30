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

// define coordinates for two-point functions
#define COOR2(x, y, size) (((x) * (size)) + (y))

/**
 * Define boundary conditions to be implemented.
 */
#define PBC
//#define APBC

/**
 * These are the scaling exponent used if SCALING_VARIABLES is (true).
 * Indeed, if SCALING_VARIABLES (true) all variables passed to the program will be scaled such that g'_i*L^{EXPONENT_g_i} is constant.
 * {g_i} -> {g'_i}={g_i_c + g_i/L^{EXPONENT_g_i}}
*/
#define SCALING_VARIABLES (true)
#define EXPONENT_G (1.0)
#define EXPONENT_LAMBDA (1.0)
#define EXPONENT_KAPPA (15.0/8.0)
#define EXPONENT_H (15.0/8.0)
#define EXPONENT_P (0.0)

/**
 * If ROUND_TRIP is defined, stop the KZ protocol when the parameters are all the same as the starting ones.
 * From here we can choose between different schemes for time evolution.
 * -  KZ_PROTOCOL : KZ of one parameter of the hamiltonian
 * - MEASUREMENT_CENTRAL_SPIN_PROTOCOL : The dynamic is triggered by monitoring continuously the central spin system.
 *                          if MEASUREMENT_CENTRAL_SPIN_PROTOCOL (true), we project the spin either on the up or down direction.
 *                          if MEASUREMENT_CENTRAL_SPIN_PROTOCOL (false), the spin is measured along the axis without projecting specifically on the up or down component
*/
//#define KZ_PROTOCOL
#define ROUND_TRIP
#define N_CYCLES (1000)
#define MEASUREMENT_PROTOCOL
// just one of these variables should be set to (true)
//****************************************************************************
// active if defind and set to (true)
#define MEASUREMENT_CENTRAL_SPIN_PROTOCOL (true)
// active if defined and set to true (true)
//#define MEASUREMENT_CHAIN_PROTOCOL (true)

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
#define MEASURE_SIGMAX
//#define MEASURE_SIGMAZ
//#define MEASURE_RANDOM


/**
 * Observables taken into account.
 * The variable WRITE_OUTPUT_EVERY_N_STEP variable is an integer defining how often the output data will be written. 
 */
#define CORRELATION_DIRECTION (1)
#define WRITE_OUTPUT_EVERY_N_STEP (10)
#define OBS_MAG
//#define OBS_ADIABATICITY
#define OBS_CORRELATION_CHAIN

/**
 * Define useful macros
 */
#define PI (3.1415926535897932384626433)
#define I (std::complex<double> (0.0, 1.0))
#define EPSILON (1e-15)

#endif
