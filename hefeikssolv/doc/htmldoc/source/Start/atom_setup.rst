.. _setup_atom:

***************
**Setup Atom**
***************

The setup for atoms in KSSOLV2.0 is simple. Here we set up H and Si atoms as
an example, which will be used in the following more complicated examples.

.. literalinclude:: ../../../../example/doc_atom_setup.m
   :language: matlab

The above example shows two ways of construction for atoms. The input for
:class:`Atom` can be either the atom number or atom symbol.

In the command window, typing a1 provides detailed information about
the atom.

.. code-block:: matlab

    >> a2

    a2 = 

      Atom with properties:
    
        symbol: 'Si'
          anum: 14
         amass: 28.0860
         venum: 4
          iloc: 1
          occs: 2
          occp: 2
          occd: 0
           iso: 0
            ic: 0
         isref: 1
         ipref: 1
         idref: 1


