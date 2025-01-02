Introduction
==================================

KSSOLV2.0 :cite:`KSSOLV2` is a MATLAB toolbox for performing electronic structure calculations
for molecules and solids at both the ground and excited states.  

************************************
**Ground state calculation**
************************************
The ground state calculation is based on the Kohn-Sham density functional 
theory.  The main equations to be solved are the Kohn-Sham equations of the form

.. math:: H(\rho) \psi_i({\bf r}) = \varepsilon_i \phi_i({\bf r}), 
   :label: kseig

where :math:`\psi_i`, :math:`i=1,2,...n_e`, with :math:`n_e` being the number 
of electrons, are orthonormal quasi-electron orbitals, and 
:math:`\varepsilon_i` are the corresponding quasi-electron energies. 
The function :math:`\rho({\bf r})` is the electron density defined to be

.. math:: \rho({\bf r}) = \sum_{i=1}^{n_e} \left| \psi_i({\bf r}) \right|^2. 

In the so called local density approximation (LDA) or general gradient approximation (GGA), the Kohn-Sham Hamiltonian :math:`H(\rho)` has the form

.. math:: H(\rho) = -\nabla^2 + \int \frac{\rho({\bf r}')}{\left| {\bf r}-{\bf r}' \right|} d{\bf r}' + v_{xc}(\rho({\bf r})) + v_{ext}({\bf r}),
   :label: ksham

where the first term in the above equation is associated with the kinectic 
energy of the electrons, the second term corresponds to the electro-static 
repulsion among electrons, the third terms is the so called 
exchange-correlation potential that accounts for other many-body properties
of the system, and the last term represents the ionic potential contributed
by nuclei. There are a number of expressions for :math:`v_{xc}` since
the exact form of the exchange-correlation potential is unknown. These 
expressions also depend on the way :math:`v_{ext}` is approximated. By 
treating inner electrons as part of an ionic core, an approach known
as the *pseudopotential* method, we obtain :math:`v_{ext}` that is easier to 
compute.

The Kohn-Sham equations form a nonlinear eigenvalue problem in which 
the Kohn-Sham Hamiltonian to be diagonalized is a function of the 
electron density :math:`\rho({\bf r})`, which is in turn a function 
of the eigenfunction :math:`\psi_i`'s to be computed.

These equations are the first order necessary condition of the 
contrained minimization problem:

.. math:: \min_{\langle\psi_i,\psi_j\rangle = \delta_{i,j}} E_{total}(\{\psi_j\}) \equiv \frac{1}{2}\sum_{i=1}^{n_e} \int \left|\nabla \psi_i ({\bf r}) \right|^2 dr + E_{ion} + E_{coul}(\rho) + E_{xc}(\rho), 
   :label: ksetot

where 

.. math:: E_{ion} = \int \rho({\bf r}) v_{ext} ({\bf r}) dr,

.. math:: E_{coul} =  \frac{1}{2}\int \int \frac{\rho(r)\rho(r')}{\left|r - r'\right|} dr dr'.

Therefore, the Kohn-Sham nonlinear eigenvalue problem can also be solved
as a constrained optimization problem.



************************************
**Excited state calculation**
************************************

The ground state is important because it determines the structure 
of a molecule or a solid (crystal strcture) at its equilibrium 
and other properties (such as elasticity) that depends on the 
structure. 
Many other properties of molecules and solids such as their ability to 
absorb photons or emit light after absorption depend on theirs excited states.

Excited states are often studied experimentally by using external 
electromagnetic radiaton (photons) to "excite" molecular or crystal
samples and measure their responses. The excited state properties can, 
in principle, be obtained from the eigenfunctions of the many-body
Schrodinger operator associted with eigenvalues above the ground state.
However, this approach is not practical due to its prohibitive computational
cost.

A more pratical approach is based on a many-body perturbation theory
in which the response of a many-body system to an external perturbing
potential is analyzed. Techniques have been developed to reduce the 
many-body problem to a single-particle equation. For ionization 
properties in which an electron is ejected, excitation energies
:math:`\varepsilon_i` can be formulated as poles of a single-particle 
Green's function :math:`G(r,r';\omega)`. They can be obtained by solving 
an eigenvalue problem of the form 

.. math:: [ H_0 + \Sigma(\varepsilon_i) ] \psi_i = \varepsilon_i \psi_i,
   :label: dyson1

where :math:`\Sigma(\omega)` is known as the self-energy, i.e., 
a potential energy operator that describes the interaction an electron 
exerts on itself through its interaction with other electrons. 

An explicit and computable form of :math:`\Sigma(\omega)` is unknown.
A widely used approximation expresses :math:`\Sigma(\omega)` as the 
convolution of the Green's function :math:`G(r,r';\omega)` to be 
determined and a screened Coulomb interaction often denoted by 
:math:`W(r,r';\omega)\equiv \int G(r,r',\omega+\omega')W(r,r',\omega') d\omega'`. Hence the name *GW approximation*.

Because both :math:`G(r,r';\omega)` and :math:`W(r,r';\omega)` depend on the
eigenvalues and eigenfunctions defined in equation :eq:`dyson1`, 
this equation should in principle be solved self-consistently.
However, this approach is very costly. An alternative that is often used 
in practice is the so called :math:`G_0 W_0` approach in which :math:`G` is
replaced by :math:`G_0`, which is constructed from the ground state Kohn-Sham
orbitals and eigenvalues directly, and the dielectric matrix used to screen
the Coulomb interaction in :math:`W` is also constructed from the ground 
states orbitals and eigenvalues.

The screened Coulomb interaction also plays an important role in describing
an optical absorption process in which an electron in an occupied orbital
absorbs incident radiation and transitions into an unoccupied orbital 
associated with a higher quasi-particle energy level as shown in the following
figure.

.. image:: figure/hvabsorb.png 

The intensity of absorption is a function of the incident radiation energy.
This function is often referred to as the optical absorption spectrum of a molecule
or a crystal as shown in the figure below. The position of each peak
corresponds to the amount of energy absorbed by an atomic system. The height of the
peak describes the likelihood of such an absorption event.

.. image:: figure/abscurves.png

The absorption energy and the likelihood of absorption (often referred to as
the oscillator strength) can be obtained by solving a non-Hermitian eigenvalue 
problem known as the Bethe-Salpeter equation.  They can also be approximated via 
time-dependent density functional theory (TDDFT).
