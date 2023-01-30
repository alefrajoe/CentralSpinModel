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
    // to avoid severe discontinuities at t=0 and trigger the dynamics we always project the state at t=0
    model.ComputeMeasurementDirection();
    model.ComputeProjectorsAlongDirection();
    #ifdef MEASUREMENT_CENTRAL_SPIN_PROTOCOL
    model.ProjectStateWithMeasurementOnCentralSpin();
    #endif
    #ifdef MEASUREMENT_CHAIN_PROTOCOL
    if(model.RandomUniformDouble() <= model.p)
    model.ProjectStateWithMeasurementsOnChain();
    #endif

    while(model.time <= model.tmax)
    {    
        // compute and write observables with a frequncy WRITE_OUTPUT_EVERY_N_STEP
        // if WRITE_OUTPUT_EVERY_N_STEP == 1 we write the observables to output after each model.deltat
        if(model.step % WRITE_OUTPUT_EVERY_N_STEP == 0)
        {
            model.ComputeObservables();
            model.WriteObservables();
        }

        // *******************   MEASUREMENT_CENTRAL_SPIN_PROTOCOL  ************************
        // if the measurement should be done, i.e.,
        // step != 0 and step % every_step_try_measurement == 0
        // note that if step == 0 and MEASUREMENT_CENTRAL_SPIN_PROTOCOL is defined model.state has already been projected
        #ifdef MEASUREMENT_CENTRAL_SPIN_PROTOCOL
        if(model.step != 0 && model.step % model.every_step_try_measurement == 0)
        if(model.RandomUniformDouble() <= model.p)
        {
            // measurement protocol
            model.ComputeMeasurementDirection();
            model.ComputeProjectorsAlongDirection();
            model.ProjectStateWithMeasurementOnCentralSpin();
        }
        #endif

        // *******************   MEASUREMENT_CHAIN_PROTOCOL  ************************
        #ifdef MEASUREMENT_CHAIN_PROTOCOL
        if(model.step != 0 && model.step % model.every_step_try_measurement == 0)
        {
            // measurement protocol
            model.ComputeMeasurementDirection();
            model.ComputeProjectorsAlongDirection();
            model.ProjectStateWithMeasurementsOnChain();
        }
        #endif

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
