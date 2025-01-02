
***********************
**Setup SCF Iteration**
***********************

The calculation of the :ref:`algo_scf` is build on the top of
:ref:`setup_molecule`. We first setup molecule and then pass the molecule to
:func:`func_scf`.

.. literalinclude:: ../../../../example/doc_scf_setup.m
   :language: matlab

Running the Matlab script provides us the following output in the command
window.

.. code-block:: matlab

    The pseudopotential for H is loaded from /Users/RyanLi/Dropbox/Research/KSSOLV/crd-kssolv/ppdata/default/H_ONCV_PBE-1.0.upf
    The pseudopotential for Si is loaded from /Users/RyanLi/Dropbox/Research/KSSOLV/crd-kssolv/ppdata/default/Si_ONCV_PBE-1.0.upf
    Warning: Renormalize starting charge 7.9452 to 8 
    > In corrcharge (line 21)
      In Ham (line 87)
      In iterinit (line 134)
      In scf (line 52)
      In doc_scf_setup (line 25) 
    Beging SCF calculation for SiH4...
    SCF iter   1:
    Rel Vtot Err    =            9.580e-02
    Total Energy    = -6.2142888183672e+00
    SCF iter   2:
    Rel Vtot Err    =            1.845e-02
    Total Energy    = -6.2284248318687e+00
    SCF iter   3:
    Rel Vtot Err    =            5.487e-03
    Total Energy    = -6.2287822079632e+00
    SCF iter   4:
    Rel Vtot Err    =            1.877e-03
    Total Energy    = -6.2288417935802e+00
    SCF iter   5:
    Rel Vtot Err    =            1.411e-04
    Total Energy    = -6.2288552002257e+00
    SCF iter   6:
    Rel Vtot Err    =            3.135e-05
    Total Energy    = -6.2288552008052e+00
    SCF iter   7:
    Rel Vtot Err    =            1.507e-06
    Total Energy    = -6.2288552013090e+00
    SCF iter   8:
    Rel Vtot Err    =            3.723e-07
    Total Energy    = -6.2288552013125e+00
    SCF iter   9:
    Rel Vtot Err    =            1.110e-08
    Total Energy    = -6.2288552013125e+00
    SCF iter  10:
    Rel Vtot Err    =            2.105e-09
    Total Energy    = -6.2288552013125e+00
    Convergence is reached!
    Etot            = -6.2288552013125e+00
    Eone-electron   = -5.3005632778731e+00
    Ehartree        =  3.2037906850856e+00
    Exc             = -2.5868682555564e+00
    Eewald          = -1.5452143529686e+00
    Ealphat         =  0.0000000000000e+00
    --------------------------------------
    Total time used =            1.139e+01
    ||HX-XD||_F     =            2.162e-09	

By default, the pseudopotential will be read from `ONCV
<http://www.quantum-simulation.org/potentials/sg15_oncv/>`_ pseudopotential.

YL: Left Here.

In this example, the ground state energy of the silane molecule which
contains one silicon atom and four hydrogen atoms is computed via a
self-consistent field (SCF) iteration. The atomic structure of the molecule
is shown below.

.. image:: figure/sih4.jpg

Let's ignore the first few lines of output for the moment. Starting from
the line ::

   Beging SCF calculation for SiH4...

the program prints out the convergence history of the self-consistent
field (SCF) iteration. In particular, it prints out the total energy
associated with the approximate solution to the Kohn-Sham equation, e.g., ::

   Total Energy    = -6.1670850083797e+00 

as well as the relative difference between the input and output potential 
in each SCF cycle ::

   Rel Vtot Err    =            8.583e-02 

defined as

:math:`\|V_{\mathrm{in}} - V_{\mathrm{out}}\|/\|V_{\mathrm{in}}\|` 

Once the SCF iteration terminates, we can examine the converged electron
density visually by typing the following at the MATLAB prompt: ::

   >> KSMolViewer(mol,H)

A GUI window will pop up. By checking the box labeled by 'density', and
dragging the bars labeled by 'Isosurface Value' and 'Transparency', we can
obtain a figure that looks like the following:

.. image:: figure/KSMolViewer.png

The molecule image can be further saved as follows.

.. image:: figure/SiH4_SCF.png

