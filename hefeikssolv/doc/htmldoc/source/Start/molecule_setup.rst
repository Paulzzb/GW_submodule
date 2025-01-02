.. _setup_molecule:

*******************
**Setup Molecule**
*******************

The setup for molecules in KSSOLV2.0 is based on :ref:`setup_atom`. Here we
set up :math:`SiH_4` molecule as an example, which will be used in the
following examples.

.. literalinclude:: ../../../../example/doc_molecule_setup.m
   :language: matlab

The above example shows a construction for molecule :math:`SiH_4`. In most
cases, the atom location must be set accurately to obtain convergence in the
:ref:`algo_scf`.

In the command window, typing mol provides detailed information about
the molecule.

.. code-block:: matlab

    >> mol

	mol = 

	  Molecule with properties:
    
               name: 'SiH4'
          supercell: [3x3 double]
              atoms: [2x1 Atom]
              alist: [5x1 double]
            xyzlist: [5x3 double]
               ecut: 12.5000
              ecut2: 50
                 n1: 32
                 n2: 32
                 n3: 32
               vext: [32x32x32 double]
              nspin: 1
        temperature: 0
             natoms: [4 1]
                vol: 1000
                nel: 8
              ppvar: []
           xyzforce: []


