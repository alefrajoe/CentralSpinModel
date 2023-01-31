# Central Spin Model

The hamiltonian of the central spin model follows

$$\hat{H}=-J\sum_{i=1}^L\hat{\sigma}^{(1)}_i\hat{\sigma}^{(1)}_{i+1}-g\sum_{i=1}^L\hat{\sigma}^{(3)}_i-h\sum_{i=1}^L\hat{\sigma}^{(1)}_i-\frac{\lambda}{2}\Sigma^{(3)}-\kappa\sum_{i=1}^L\hat{\Sigma}^{(a)}\hat{\sigma}_i^{(b)}$$

where Periodic Boundary Conditions (PBC) are intended in the above expression, and $a, b=1,2,3$ are two directions chosen at running time. We set $J=1$, so that all quantities are measured in units of $J$.

To compile the code

> g++ -O3 *.cpp -o central -lm -larmadillo -lstdc++fs

In *macro.hpp*, if the variable SCALING_VARIABLES is set to **true**, all those parameters that in *macro.hpp* have the associated
variable #define EXPONENT_.. (double) defined are rescaled with the same exponent contained within the parentheses. Example

> #define SCALING_VARIABLES (true)
> #define EXPONENT_H (15.0/8.0)

The value of h passed is rescaled by a factor $h \to h^\prime=h/L^{15/8}$.

## Measurement induced dynamic

To use a measurement protocol, the proper macro should be defined in *macro.hpp*.

> #define MEASUREMENT_PROTOCOL

Two distinguished measurement-induced-dynamic protocols can be devised with this code.
On the one hand, the central spin qubit can be continuously monitored, projecting its state along the up or down component along direction $\vec{n}$ (the direction of the observation). In this case, the spin is observed after a time frame $t_M$ with a probability $p$. To use this protocol set

> #define MEASUREMENT_CENTRAL_SPIN_PROTOCOL (true)

On the other hand, the Ising chain can be continuously monitored, given a probability $p$ per site of observation. In this case, instead, each one among the $L$ spins has a probability $p$ of being projecting along the up or down component of the direction observed. See [RV-2020](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.102.035119).

> #define MEASUREMENT_CHAIN_PROTOCOL (true)

Three possible spin projection are feasible within these two protocols.
To measure the spin(s) along the $x$ axis.

> #define MEASURE_SIGMAX

To measure the spin(s) along the $z$ axis.

> #define MEASURE_SIGMAZ

To measure the spin(s) along a random axis.

> #define MEASURE_RANDOM

The parameters to be passed to the program are
- **L** (**L**): int - length of the Ising chain
- **g** (**g**): double - transverse coupling of the chain
- **h** (**h**): double - longitudinal coupling of the chain
- **$\lambda$** (**lambda**): double - gap of the single qubit
- **$\kappa$** (**kappa**): double - interaction strength between the qubit and the chain
- **a** (**a**): int - interaction term for the qubit $\sim \hat{\Sigma}^{(a)}$, it can be $a=1,2,3$.
- **b** (**b**): int - interaction term for the chain $\sim \hat{\sigma}^{(b)}$, it can be $b=1,2,3$.
- **t_M** (**tm**): double - time frame after which a measurement is tried.
- **p** (**p**): double - probability of projecting the spin(s), after a time frame $t_M$, along the direction chosen in *macro.hpp*.
- **seed** (**seed**): int - seed used for the random number generator. 
- **steptm** (**steptm**): int - the integer is used to determing **$\delta t$** used in the simulation. If passed to the program, $\delta t=t_M / steptm$. 
- **$t_{MAX}$** (**tmax**): double - diration of the time evolution.

## Kibble-Zurek protocol

The parameters to be passed to the program are
- **L** (**L**): int - length of the Ising chain
- **g** (**g**): double - transverse coupling of the chain
- **h** (**h**): double - longitudinal coupling of the chain
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

After compiling the code as described above, you can run a simulation passing the parameter name (the one inside the round brackets) and the value for that observables. For instance, to run the program with $L=6, g=1.1, \lambda=0.5, \kappa=1, a=1, b=1, t_{KZ}=50, \delta t=0.002, end=1.5$

> ./central L 6 g 1.1 lambda 0.5 kappa 1 a 1 b 1 tkz 50 dt 0.002 end 1.5

# Observables
We focus on the density matrix of the central spin to study our system.
To this purpose, we first point out that the density matrix $\rho$ can be written as follows
$$\rho=\frac{1}{2}\bigg(\mathbb{1}+\langle\vec{\Sigma}\rangle\cdot\vec{\Sigma}\bigg)$$.
Thus all quantities for the central spin can be obtained from the magnetizations along the three directions $\langle\vec{\Sigma}\rangle$.
The simulation returns
- **Central spin magnetization** $\langle \Sigma^{(i)}(t)\rangle$:
> \# magx, magy, magz

$$\langle \Sigma^{(i)} \rangle, \ \text{for} \ i=1,2,3$$.

Magnetizations along the three directions of the current state.
- **Ground state Central spin magnetization** $\langle \Sigma^{(i)}_{GS}(t)\rangle$:
> \# magGSx, magGSy, magGSz

$$\langle \Sigma^{(i)}_{GS} \rangle, \ \text{for} \ i=1,2,3$$

Magnetizations along the three directions for the ground state corresponding to the current hamiltonian paramters.
- **Adiabaticity** $A(t)$:
> \# adiabaticity

$$A(t) = |\braket{\Omega_{\{C_i\}}|\psi(t)}|$$
where $\ket{\Omega_{\{C_i\}}}$ is the ground state of the hamiltonian with parameters $\{C_i\}$, i.e., the current hamiltonian.
- **Decoherence** $D(t)$
$$D(t) = 1 - \text{Tr}\rho^2(t)$$
where $\rho$ is the density matrix of the central spin. The decoherence is not computed directly by this code.
However, note that
$$D(t)=\frac{1}{2}(1-\langle\vec{\Sigma}\rangle^2)$$

- **Longitudinal Ising chain magnetization** $\frac{1}{L}\displaystyle\sum_{x=1}^L \langle \sigma^{(i)}_{x}\rangle$:
> \# magChain

$$\frac{1}{L}\sum_{x=1}^L\langle \sigma^{(i)}_{x} \rangle, \ \text{for} \ i=1$$

- **Longitudinal Ising chain (disconnected) susceptibility** :
> \# suscChain

$$\frac{1}{L}\sum_{x, y=1}^L\langle \sigma^{(i)}_{x}\sigma^{(i)}_{y} \rangle, \ \text{for} \ i=1$$