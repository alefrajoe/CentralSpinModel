#ifndef HAMILTONIAN_HPP
#define HAMILTONIAN_HPP

#include "macro.hpp"

class Model
{
    public:
    // hamiltonian parameters
    int L;
    double g;
    double lambda;
    double kappa;
    double h;
    double deltat;
    double time;
    double t_KZ;
    double final_param;
    double p;
    double tm;
    double tmax;
    int step;
    int every_step_try_measurement;
    int interaction_spin;
    int interaction_chain;

    // random number generator
    std::mt19937 random_engine;
    std::uniform_real_distribution<double> uniform_distribution{0.0, 1.0};
    int seed;
    // random measurements
    double nx;
    double ny;
    double nz;
    arma::sp_cx_dmat id22;
    arma::sp_cx_dmat up_proj;
    arma::sp_cx_dmat down_proj;

    // operators related to the magnetization of the central spin 
    arma::sp_cx_dmat *magx;
    arma::sp_cx_dmat *magy;
    arma::sp_cx_dmat *magz;

    // groundstate properties
    double eigenvalues[2];
    double adiabaticity;
    arma::cx_vec *state;

    // observables
    double magObs[3];
    double maggroundObs[3];

    // filename and hamiltonian
    std::string filename;
    arma::sp_cx_dmat *hamiltonian;
    std::string real_or_complex;

    // initializer
    Model(int argc, char **argv);

    // hamiltonian terms
    void AddSpinHamiltonian(double par);
    void AddTransverseFieldChain(double par);
    void AddLongitudinalFieldChain(double par);
    void AddLongitudinalHoppingChain();
    void AddInteractionCentralSpinAndChain(int a, int b, double par);
    void AddHamiltonian();

    // diagonalization
    void GroundStateAndEigenvals(arma::cx_vec *vec, bool replace_eigvals);
    // state normalization
    void StateNormalization();

    // random numbers
    double RandomUniformDouble();

    // compute or write observables
    double ExpectationValueOfOperatorOnState(arma::sp_cx_dmat *op, arma::cx_vec *vec);
    double ComputeAdiabaticity();
    void ComputeObservables();
    void WriteObservables();

    // time evolution
    void EvolveHamiltonianByDx(double dx);
    void RungeKuttaStep();
    void TimeEvolutionProtocol();

    // measurement induced dynamics
    void ComputeMeasurementDirection();
    void ComputeProjectorsAlongDirection();
    void ProjectStateWithMeasurementOnCentralSpin();

    // projectors along the chain
    int proj_chain_count;
    std::vector<arma::sp_cx_dmat> *proj_chain;
    void InitializeProjChain();
    void SingleStepConstructIterativeProjectorsChain();
    void ProjectStateWithMeasurementsOnChain();
};

#endif