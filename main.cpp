#include "macro.hpp"
#include "hamiltonian.hpp"

using namespace std;
using namespace arma;

int main(int argc, char **argv)
{
    std::cout << std::setprecision(16); 

    // *********** Read parameters from input **********
    // read parameters from **argv
    // create the model
    Model model(argc, argv);

    // ********** construct the hamiltonian and compute the Ground State ************* 
    model.AddHamiltonian();
    model.GroundStateAndEigenvals(model.state);

    // ************   Kibble - Zurek  **************************
    model.TimeEvolutionProtocol();

    // return exit success
    return 0;
}
