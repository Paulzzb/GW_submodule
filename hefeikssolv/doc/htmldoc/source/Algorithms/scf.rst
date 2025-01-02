.. _algo_scf:

Self Consistent Field Iteration
===============================

After the Kohn-Sham equations are discretized by planewaves,
they become a set of nonlinear algebraic
eigenvalue equations 

.. math:: H(\rho(X)) X = X\Lambda  

in which the matrix to be diagonalized
is a function of the eigenvectors :math:`X` to be computed.
The desire eigenvalues to be computed (on the diagonal of :math:`\Lambda`) 
are the :math:`n_e` smallest eigenvalues of :math:`H`, where :math:`n_e`
is proportional to the number of electrons.
Due to this nonlinear coupling between the matrix operator and
its eigenvectors, the Kohn-Sham equations are much more difficult to
solve than a standard linear eigenvalue problem.

Currently, the most widely used numerical method for solving the
this type of problem is the Self Consistent Field (SCF) iteration.
An alternative approach which is designed to minimize the
Kohn-Sham total energy directly will be reviewed later.
Both methods have been implemented in KSSOLV to allow users to
compute the solution to the Kohn-Sham equations associated with various
molecules and solids.

To use the SCF algorithm in KSSOLV, we can simply type::

  >> [mol, info] = scf(mol);

where mol is a previously constructed *Molecule* object.
The function call returns a *Molecule* object which may contain updated
atomic forces.  The total energy of the molecule obtained at the end of the 
SCF procedure, as well as the eigenvalues and eigenvectors of the Kohn-Sham
Hamiltonian are stored in the info structure.

Given an initial guess of the eigenvector matrix :math:`X` and electron 
density :math:`\rho`, the simplest SCF procedure can be described by the 
following simple loop ::
  
  iter = 1
  while (not converged and iter < maxiter)
     Compute the ne smallest eigenvalues and corresponding eigenvector of H
     Update the eletron density, potential and H
     iter = iter + 1
  end while

However, this simple procedure rarely converges because it corresponds to
a simple fixed point iteration for a set of nonlinear equations satisfied by
the electron density

.. math::  \rho = f_\mu (\rho), \sum_i \rho(i) = n_e,

where the function :math:`f_\mu(\rho)` can be written as the diagonal part of 
a step function applied to :math:`H(\rho)`. The step function assume the 
value of 1 on the interval :math:`(-\infty,\mu)` and 0 on the interval 
:math:`[\mu,\infty]`. The value of :math:`\mu` is called the chemical potential. It is not uniquely defined by the above equation if there is a gap between the :math:`n_e`-th and the :math:`n_e+1`-st eigenvalues of the converged :math:`H(\rho)`.
          
Because :math:`f_\mu(\rho)` is generally not a global contraction, there is
no guarantee that a simple fixed point iteration that starts from any initial
guess of :math:`\rho` can converge. A more sophisticated iterative algorithm
is required to solve the nonlinear system successfully. We will discuss this 
alogrithm in the next section.


