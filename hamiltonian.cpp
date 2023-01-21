#include "hamiltonian.hpp"

/**
 * Model constructor.
 * ------------------------------------
 * parameters:
 *              - int argc
 *              - char **argv
 * return:
 *              - None
*/
Model::Model(int argc, char **argv)
{
    // initialize all parameters required for the simulation
    this->L = atoi(argv[1]);
    this->g = atof(argv[2]);
    this->lambda = atof(argv[3]);
    this->kappa = atof(argv[4]);
    this->interaction_spin = atoi(argv[5]);
    this->interaction_chain = atoi(argv[6]);
    
    // initialize the hamiltonian
    // the 1^st spin is the central qubit
    this->hamiltonian = new arma::sp_cx_dmat(pow(2, this->L+1), pow(2, this->L+1));

    // initialize the mag* operators ----------------------------------------
    this->magx = new arma::sp_cx_dmat(pow(2, this->L+1), pow(2, this->L+1));
    this->magy = new arma::sp_cx_dmat(pow(2, this->L+1), pow(2, this->L+1));
    this->magz = new arma::sp_cx_dmat(pow(2, this->L+1), pow(2, this->L+1));

    // define Pauli matrices
    arma::sp_cx_dmat sigmax(2, 2), sigmay(2, 2), sigmaz(2, 2), id(pow(2, this->L), pow(2, this->L));
    sigmax(0, 1) = 1.0;     // sigmax
    sigmax(1, 0) = 1.0;
    sigmay(0, 1) = -1.0j;       // sigmay
    sigmay(1, 0) = 1.0j;
    sigmaz(0, 0) = 1.0;         // sigmaz
    sigmaz(1, 1) = -1.0;

    // construct mag* operators
    (*magx) = arma::kron(sigmax, id.eye(pow(2, this->L), pow(2, this->L)));
    (*magy) = arma::kron(sigmay, id.eye(pow(2, this->L), pow(2, this->L)));
    (*magz) = arma::kron(sigmaz, id.eye(pow(2, this->L), pow(2, this->L)));

    // initialize state to nothing
    this->state = new arma::cx_vec(pow(2, this->L+1));
    
    // initialize observables
    for(int i=0; i<3; i++) this->magObs[i] = 0.0;

    // create the directory where data will be saved (if it doesn't exist)
    std::string directory{"data_centralspin"};
    std::experimental::filesystem::create_directory(directory);

    // create substring for the interaction
    std::string inter_string{""};
    switch (this->interaction_spin)
    {
    case 1:
        inter_string+="X";
        break;
    case 2:
        inter_string+="Y";
        break;
    case 3:
        inter_string+="Z";
        break;
    default:
        std::cout << "Error interaction spin passed!" << std::endl;
        exit(1);
    }
    switch (this->interaction_chain)
    {
    case 1:
        inter_string+="X";
        break;
    case 2:
        inter_string+="Y";
        break;
    case 3:
        inter_string+="Z";
        break;
    default:
        std::cout << "Error interaction chain passed!" << std::endl;
        exit(1);
    }
    // create filename
    this->filename = "data_centralspin/data" + inter_string + std::to_string(this->g) + "g" + std::to_string(this->lambda) + "lambda" + std::to_string(this->kappa) + "kappa" + std::to_string(this->L) + "L.txt";

    // open file inside directory and write first line
    // create ofstream variable
    std::ofstream outfile;
    // open file
    outfile.open(this->filename, std::ios_base::app); 
    outfile << "#   L   g   lambda  kappa   magx    magy    magz" << std::endl;
    // close file
    outfile.close();
}

/**
 * Add the hamiltonian of the central spin H_q to this->hamiltonian.
 * ------------------------------------
 * parameters:
 *              - double : par  - the hamiltonian term is multiplied by (par / 2.0)
*/
void Model::AddSpinHamiltonian(double par)
{
    // define sigmaz
    arma::sp_cx_dmat sigmaz(2, 2);
    sigmaz(0, 0) = 1.0;
    sigmaz(1, 1) = -1.0;

    // define identity
    arma::sp_cx_dmat id(pow(2, this->L), pow(2, this->L));
    id.eye(pow(2, this->L), pow(2, this->L));

    // add the term to hamiltonian
    (*this->hamiltonian) = (*this->hamiltonian) - (par / 2.0) * arma::kron(sigmaz, id);
}

/**
 * Add the transverse field hamiltonian to the Ising chain to this->hamiltonian.
 * ------------------------------------
 * parameters:
 *              - double : par  - the hamiltonian term is multiplied by par
*/
void Model::AddTransverseFieldChain(double par)
{
    // initialize 2 X 2 identity and sigmaz
    arma::sp_cx_dmat id(2, 2);
    id.eye(2, 2);
    arma::sp_cx_dmat sigmaz(2, 2);
    sigmaz(0, 0) = 1.0;
    sigmaz(1, 1) = -1.0;

    // add all L terms to the Ising chain
    for(int i=0; i<this->L; i++)
    {
        // initialize temp to 2 X 2 identity
        arma::sp_cx_dmat temp(2, 2);
        temp.eye(2, 2);

        // generate the hamiltonian term
        for(int j=0; j<this->L; j++)
        {
            // if i == j kron with sigmax
            if(i == j) temp = arma::kron(temp, sigmaz);
            // else kron with id
            else temp = arma::kron(temp, id);
        }

        // add the temp term to the hamiltonian
        (*this->hamiltonian) = (*this->hamiltonian) - ((par) * temp);
    }
}

/**
 * Add the longitudinal hopping term for the Ising chain to this->hamiltonian.
 * ------------------------------------
 * parameters:
 *              - None
 * return:
 *              - None
*/
void Model::AddLongitudinalHoppingChain()
{
    // initialize 2 X 2 identity and sigmaz
    arma::sp_cx_dmat id(2, 2);
    id.eye(2, 2);
    arma::sp_cx_dmat sigmax(2, 2);
    sigmax(0, 1) = 1.0;
    sigmax(1, 0) = 1.0;

    // add all L-1 terms to the Ising chain
    for(int i=0; i<this->L-1; i++)
    {
        // initialize temp to 2 X 2 identity (this operator acts on the central spin)
        arma::sp_cx_dmat temp(2, 2);
        temp.eye(2, 2);

        // generate the hamiltonian term
        for(int j=0; j<this->L; j++)
        {
            // if j == i or j == i+1 kron with sigmax
            if((j == i) || (j == i+1)) temp = arma::kron(temp, sigmax);

            // else kron with id
            else temp = arma::kron(temp, id);
        }

        // add the temp term to the hamiltonian
        (*this->hamiltonian) = (*this->hamiltonian) - (temp);
    }    

    // add the last term depending on the boundary conditions
    // initialize temp to 2 X 2 identity (this operator acts on the central spin)
    arma::sp_cx_dmat temp(2, 2);
    temp.eye(2, 2);
    for(int i=0; i<this->L; i++)
    {
        // add the kinetic term between site 0 and L-1
        if(i == 0 || i == this->L-1) temp = arma::kron(temp, sigmax);
        // else multiply by the identity
        else temp = arma::kron(temp, id);
    }

    #ifdef PBC
    (*this->hamiltonian) = (*this->hamiltonian) - temp;
    #endif
    #ifdef APBC
    (*this->hamiltonian) = (*this->hamiltonian) + temp;
    #endif
}

/**
 * Add the interaction term between the Ising chain and the central spin.
 * The interaction term is of the form \sum_{i=1}^L \Sigma^(a)\sigma_i^(b), where a, b = 1, 2, 3.
 * ------------------------------------
 * parameters:
 *              - int a : direction along which the central spin interaction is oriented
 *              - int b : direction along whihc the chain interaction is oriented
 *              - double par : coupling strength of the interaction considered
*/
void Model::AddInteractionCentralSpinAndChain(int a, int b, double par)
{
    // define 2 X 2 identity matrix
    arma::sp_cx_dmat id(2, 2);
    id.eye(2, 2);

    // first, define the three Pauli matrices
    arma::sp_cx_dmat sigmax(2, 2), sigmay(2, 2), sigmaz(2, 2);
    // sigmax
    sigmax(0, 1) = 1.0;
    sigmax(1, 0) = 1.0;
    // sigmay
    sigmay(0, 1) = -1.0j;
    sigmay(1, 0) = 1.0j;
    // sigmaz
    sigmaz(0, 0) = 1.0;
    sigmaz(1, 1) = -1.0;

    // initialize interaction variables
    arma::sp_cx_dmat interaction_spin, interaction_chain;
    // assign interaction term for the central spin
    switch (a)
    {
        case 1:
            interaction_spin = sigmax;
            break;
        case 2:
            interaction_spin = sigmay;
            break;
        case 3:
            interaction_spin = sigmaz;
            break;
        default:
            exit(1);
            break;
    }
    // assign interaction term for the chain
    switch (b)
    {
        case 1:
            interaction_chain = sigmax;
            break;
        case 2:
            interaction_chain = sigmay;
            break;
        case 3:
            interaction_chain = sigmaz;
            break;
        default:
            exit(1);
            break;
    }

    // construct interaction
    // for L terms in the Ising chain
    for(int i=0; i<this->L; i++)
    {
        // initialize temp variable
        arma::sp_cx_dmat temp(2, 2);
        temp = (interaction_spin);

        // cycle over L terms in the Ising chain
        for(int j=0; j<this->L; j++)
        {
            // if i == j
            if(i == j) temp = arma::kron(temp, interaction_chain);
            // else id
            else temp = arma::kron(temp, id);
        }

        // add this term to the hamiltonian
        (*this->hamiltonian) = (*this->hamiltonian) - (par * temp);
    }
}

/**
 * Generate the hamiltonian of the model.
 * No arguments are required since they are all contained into the model.
 * ------------------------------------
 * parameters:
 *              - None
 * return:
 *              - None
*/
void Model::AddHamiltonian()
{
    // add all terms of the hamiltonian
    this->AddSpinHamiltonian(this->lambda);
    this->AddTransverseFieldChain(this->g);
    this->AddLongitudinalHoppingChain();
    this->AddInteractionCentralSpinAndChain(this->interaction_spin, this->interaction_chain, this->kappa);
}

/**
 * Compute the groundstate of the Hermitian hamiltonian and the first 2 eigenvals (these are sufficient to compute the gap).
 * ------------------------------------
 * parameters:
 *              - arma::cx_vec *vec : vector where the groundstate of the hamiltonian will be saved
*/
void Model::GroundStateAndEigenvals(arma::cx_vec *vec)
{
    // auxiliary matrix
    arma::cx_vec eigenvalues;
    arma::cx_mat eigvec;

    // compute lowest eigenvalues and eigenvectors
    arma::eigs_gen(eigenvalues, eigvec, (*this->hamiltonian), 2, "sr");

    // save the groundstate
    (*vec) = eigvec.col(0);

    // save eigenvalues into the model
    this->eigenvalues[0] = eigenvalues.at(0).real();
    this->eigenvalues[1] = eigenvalues.at(1).real();

    // print gap
    std::cout << "gap is " << this->eigenvalues[1] - this->eigenvalues[0] << std::endl; 
}

/**
 * Compute the expectation value of the operator passed to the method on the state "this->state".
 * The function returns a real number.
 * ------------------------------------
 * parameters:
 *              - arma::sp_cx_dmat *op : operator whose expectation value will be evaluated < *op >
 * return:
 *              - double : return a real number equal to <state| *op | state> 
*/
double Model::ExpectationValueOfOperatorOnState(arma::sp_cx_dmat *op)
{
    std::complex<double> expval{0.0};
    
    // compute temp = *op | state> 
    arma::cx_vec temp = (*op) * (*this->state);

    // compute <state | temp>
    for(int i=0; i<pow(2, this->L+1); i++) expval += conj(this->state->at(i)) * temp.at(i);

    return expval.real();
}

/**
 * Compute all observables on the current state.
 * The observables to be observed are contrelled from the macro file.
 * ------------------------------------
 * parameters:
 *              - None
 * return:
 *              - None
*/
void Model::ComputeObservables()
{
    #ifdef OBS_MAG
    this->magObs[0] = this->ExpectationValueOfOperatorOnState(this->magx);
    this->magObs[1] = this->ExpectationValueOfOperatorOnState(this->magy);
    this->magObs[2] = this->ExpectationValueOfOperatorOnState(this->magz);
    #endif
}

/**
 * Write all the observables in the txt file.
 * The file is closed before the function ends.
 * ------------------------------------
 * parameters:
 *              - None
 * return:
 *              - None
*/
void Model::WriteObservables()
{
    // create ofstream variable
    std::ofstream outfile;
    
    // open file
    outfile.open(this->filename, std::ios_base::app); 

    // close file
    outfile.close();
}