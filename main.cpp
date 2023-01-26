#include "macro.hpp"
#include "hamiltonian.hpp"

using namespace std;
using namespace arma;

int main(int argc, char **argv)
{
    // set stdout precision to 16 digits
    std::cout << std::setprecision(16);

    // *********** Read parameters from input **********
    // read parameters from **argv
    // create the model
    Model model(argc, argv);

    // ********** construct the hamiltonian and compute the Ground State ************* 
    model.AddHamiltonian();
    model.GroundStateAndEigenvals(model.state, true);

    // ************   Kibble - Zurek  **************************
    for(int i=0; i<N_CYCLES; i++) model.TimeEvolutionProtocol();
    // compute observables also at the end of the simulation
    model.ComputeObservables();
    // write to file also at the end of the simulation
    model.WriteObservables();

    // return exit success
    return 0;
}
