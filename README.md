# Central Spin Model
The hamiltonian of the central spin model follows

$$\hat{H}=-J\sum_{i=1}^L\hat{\sigma}^{(1)}_i\hat{\sigma}^{(1)}_{i+1}-g\sum_{i=1}^L\hat{\sigma}^{(3)}_i-\frac{\lambda}{2}\Sigma^{(3)}-\kappa\sum_{i=1}^L\hat{\Sigma}^{(a)}\hat{\sigma}_i^{(b)}$$

where Periodic Boundary Conditions (PBC) are intended, and $a, b=1,2,3$ are two directions chosen at running time.
We set $J=1$, so that all quantities are measured in units of $J$.

To compile the code
> g++ -O3 *.cpp -o central -lm -larmadillo -lstdc++fs

The parameters to be passed to the program are
- **L** (**L**): int - length of the Ising chain
- **g** (**g**): double - transverse field coupling of the chain
- **$\lambda$** (**lambda**): double - gap of the single qubit
- **$\kappa$** (**kappa**): double - interaction strength between the qubit and the chain
- **a** (**a**): int - interaction term for the qubit $\sim \hat{\Sigma}^{(a)}$, it can be $a=1,2,3$.
- **b** (**b**): int - interaction term for the chain $\sim \hat{\sigma}^{(b)}$, it can be $b=1,2,3$.
- **$t_{KZ}$** (**tkz**): double - typical timescale over which the coupling in the KZ protocol is changed.
                        The coupling could be $C=g,\lambda$ or $\kappa$, and it is changed linearly with time as
                        $C\to C+t/t_{KZ}$.
- **$\delta$ t** (**dt**): double - infinitesimal time step for the $4^{th}$-order the Runge-Kutta algorithm. It should always
                                    be a positive number independently of the parameter end.
- **end** (**end**): double - last coupling to be simulated. The time evolution stops as soon as the coupling $C$ is equal to 
                    this value. The coupling changed during the KZ protocol is controlled in *macro.h*.
                    The value of the end parameter could be either larger or smaller than the starting parameter.

After compiling the code as described above, you can run a simulation passing the parameter name (the one inside the round brackets) and the value for that observables. For instance, to run the program with $L=6, g=1.1, \lambda=0.5, \kappa=1, a=1, b=1, t_{KZ}=50, \delta t=0.002, end=1.5$ launch

> ./central L 6 g 1.1 lambda 0.5 kappa 1 a 1 b 1 tkz 50 dt 0.002 end 1.5

# Observables
We consider the following observables
- **Central spin magnetization** $m^{(i)}(t)$:
$$\langle \Sigma^{(i)} \rangle, \ \text{for} \ i=1,2,3$$
- **Adiabaticity** $A(t)$:
$$A(t) = $\abs{\bra{\Omega_{\{C_i\}}}\ket{\psi(t)}}$
where $\ket{\Omega_{\{C_i\}}}$ is the ground state of the hamiltonian with parameters $\{C_i\}$, i.e., the current hamiltonian.
- **Decoherence** $D(t)$
$$D(t) = 1 - \text{Tr}\rho^2(t)$$
where $\rho$ is the density matrix of the central spin.
