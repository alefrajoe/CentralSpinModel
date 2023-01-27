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

    #ifdef MEASUREMENT_PROTOCOL
    while(model.time <= model.tmax)
    {    
        // compute and write observables
        model.ComputeObservables();
        model.WriteObservables();

        // if the measurement should be done, i.e.,
        // it is not the first step and step % every_step_try_measurement == 0
        if(model.step != 0 && model.step % model.every_step_try_measurement == 0)
        if(model.RandomUniformDouble() <= model.p)
        {
            // measurement protocol
            model.ComputeMeasurementDirection();
            model.ComputeProjectorsAlongDirection();
            model.ProjectStateWithMeasurement();
        }

        // evolve in time
        // ! StateNormalization is inside the time evolution function
        model.RungeKuttaStep();
    }
    #endif

    #ifdef KZ_PROTOCOL
    // ************   Kibble - Zurek  **************************
    for(int i=0; i<N_CYCLES; i++) model.TimeEvolutionProtocol();
    // compute observables also at the end of the simulation
    model.ComputeObservables();
    // write to file also at the end of the simulation
    model.WriteObservables();
    #endif

    // return exit success
    return 0;
}
