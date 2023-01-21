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
    double deltat;
    int interaction_spin;
    int interaction_chain;

    // operators related to the magnetization of the central spin 
    arma::sp_cx_dmat *magx;
    arma::sp_cx_dmat *magy;
    arma::sp_cx_dmat *magz;

    // groundstate properties
    double eigenvalues[2];
    arma::cx_vec *state;

    // observables
    double magObs[3];

    // filename and hamiltonian
    std::string filename;
    arma::sp_cx_dmat *hamiltonian;

    // initializer
    Model(int argc, char **argv);

    // hamiltonian terms
    void AddSpinHamiltonian(double par);
    void AddTransverseFieldChain(double par);
    void AddLongitudinalHoppingChain();
    void AddInteractionCentralSpinAndChain(int a, int b, double par);
    void AddHamiltonian();

    // diagonalization
    void GroundStateAndEigenvals(arma::cx_vec *vec);

    // compute or write observables
    double ExpectationValueOfOperatorOnState(arma::sp_cx_dmat *op);
    void ComputeObservables();
    void WriteObservables();
};

#endif