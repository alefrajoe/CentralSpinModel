# Central Spin Model
The hamiltonian of the central spin model follows

$$\hat{H}=-J\sum_{i=1}^L\hat{\sigma}^{(1)}_i\hat{\sigma}^{(1)}_{i+1}-g\sum_{i=1}^L\hat{\sigma}^{(3)}_i-\frac{\lambda}{2}\Sigma^{(3)}-\kappa\sum_{i=1}^L\hat{\Sigma}^{(a)}\hat{\sigma}_i^{(b)}$$

where Periodic Boundary Conditions (PBC) are intended, and $a, b=1,2,3$ are two directions chosen at running time.
We set $J=1$, so that all quantities are measured in units of $J$.

To compile the code
> g++ -O3 *.cpp -o central -lm -larmadillo -lstdc++fs

The parameters to be passed to the program are
- **L**: int - length of the Ising chain
- **g**: double - transverse field coupling of the chain
- **$\lambda$**: double - gap of the single qubit
- **$\kappa**: double - interaction strength between the qubit and the chain
- **a**: int - interaction term for the qubit $\sim \hat{\Sigma}^{(a)}$, it can be $a=1,2,3$.
- **b**: int - interaction term for the chain $\sim \hat{\sigma}^{(b)}$, it can be $b=1,2,3$.
- **$t_{KZ}$**: double - typical timescale over which the coupling in the KZ protocol is changed.
                        The coupling is $C=g,\lambda, \text{or} \kappa$, and it is changed with time as
                        $C\to C+t/t_{KZ}$. 

# Observables
We consider the following observables
- **Central spin magnetization** $m^{(i)}(t)$:
$$\langle \Sigma^{(i)} \rangle, \ \text{for} \ i=1,2,3$$
- **Decoherence** $D(t)$
$$D(t) = 1 - \text{Tr}\rho^2(t)$$
where $\rho$ is the density matrix of the central spin.
