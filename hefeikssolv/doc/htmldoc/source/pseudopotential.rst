Pseudo Potential
==================================

The concept of the pseudo potential is to approximate the effects of the all
electrons by a approximated system with valence electrons only.
Mathematically, we could view this idea as a projection from a high
dimensional space (the space corresponding to all electrons) down to a low
dimensional one (the space corresponding to valence electrons). And the lost
of the projection is marginal. However, the Kohn-Sham equation governing the
system must be modified as well. 

.. toctree::
   :maxdepth: 1

   pseudopot/localpot
   pseudopot/nonlocalpot


Pseudo Potential Transformation
************************************

Let us conceptually derive the pseudo potential transformation without dive
deep into the quantum calculation. We simply denote the original Hamiltonian
as :math:`\mathcal{H}`, the core states as :math:`\{\Ket{\chi_n}\}_n` and
the core eigenvalues as :math:`\{E_n\}_n`. Without loss of generality, we
assume :math:`\Ket{\psi}` is a single valence state, and the expansion on
the core states and a smoother "state" :math:`\Ket{\phi}` is given as

.. math:: \Ket{\psi} = \Ket{\phi} + \sum_n a_n\Ket{\chi_n}.
   :label: eq:valence-exp

Due to the fact that the valence state is orthogonal to all core states, we
could derive all the coefficients in :eq:`eq:valence-exp` as follows,

.. math::
   :label: eq:valence-innerprod

    \BraKet{\chi_m}{\psi} = \BraKet{\chi_m}{\phi}
    + \sum_n a_n \BraKet{\chi_m}{\chi_n} = 0.

Based on :eq:`eq:valence-innerprod`, the explicit expression for each
:math:`a_n` can be derived. And equation :eq:`eq:valence-exp` is rewritten
as,

.. math::
   :label: eq:valence-explicit

   \Ket{\psi} = \Ket{\phi} + \sum_n \BraKet{\chi_n}{\phi}\Ket{\chi_n}.

Applying the Hamiltonian :math:`\mathcal{H}` to the equation
:eq:`eq:valence-explicit`, and reorganizing the equation, we obtain,

.. math::
   :label: eq:hamiltonian-expansion

   \mathcal{H}\Ket{\phi} + \sum_n (E-E_n) \Ket{\chi_n} \BraKet{\chi_n}{\phi}
   = E \Ket{\phi}.

Equation :eq:`eq:hamiltonian-expansion` is equivalent to a updated
Hamiltonian :math:`\mathcal{H}+V=\mathcal{H}
+ \sum_n (E-E_n)\KetBra{\chi_n}{\chi_n}` operating on the smooth part of the
valence state, i.e., the :math:`\Ket{\phi}` of the :math:`\Ket{\psi}`.
And the updated system is now,

.. math::
   
   \left(\mathcal{H}+V\right)\Ket{\phi} = E \Ket{\phi}.

In this new system, the eigenstates are much smoother than the original
valence states.

Norm-Conserving Pseudo Potentials
************************************
Key properties of the pseudo potential are the optimum smoothness and
the transferability. However, these properties require four configurations:

1. All-electron and pseudo eigenvalues agree for the reference
   configuration.

   .. math::
      \mathcal{H}\Ket{\psi^{AE}} =& \epsilon \Ket{\psi^{AE}}\\
      \left(\mathcal{H}+V\right)\Ket{\psi^{PS}} =& \epsilon \Ket{\psi^{PS}}
   
   Here :math:`\Ket{\psi^{AE}}` denotes the all-electron eigenstate and
   :math:`\Ket{\psi^{PS}}` denotes the smooth eigenstate.

2. AE and PS wave functions agree beyond a certain cutoff, :math:`r_c`.

   .. math::
      \psi^{AE}(r) = \psi^{PS}(r)\quad\text{for}\quad r\geq r_c.

3. Real and pseudo norm squares integrated from :math:`0` to :math:`R` for
   all :math:`R\le r_c` agree,

   .. math::
      \int_0^R \lvert \phi^{AE}\rvert^2 r^2 \mathrm{d} r = \int_0^R \lvert
      \phi^{PS}\rvert^2 r^2 \mathrm{d} r.

4. The logarithmic derivatives and the energy derivative of the logarithmic
   derivative agree for all :math:`R\le r_c`,

   .. math::
      \left[\left(r\phi^{AE}(r)\right)^2\frac{\dif}{\dif E}
      \frac{\dif}{\dif r}\ln{\phi^{AE}(r)}\right]_R =
      \left[\left(r\phi^{PS}(r)\right)^2\frac{\dif}{\dif E}
      \frac{\dif}{\dif r}\ln{\phi^{PS}(r)}\right]_R.




