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
    model.AddHamiltonian();
    model.GroundStateAndEigenvals(model.state);
    model.ComputeObservables();
    model.WriteObservables();

    // return exit success
    return 0;
}
