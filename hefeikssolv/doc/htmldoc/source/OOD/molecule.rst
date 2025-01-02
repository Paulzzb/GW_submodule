
Molecule
==================================

Molecule is a class that defines a number of attributes of a molecule. 

Class Function List
----------------------

.. class:: Molecule

   .. method:: Molecule(str1,field1,str2,field2,...)

      Constructor for Molecule based on the field values corresponding to
      the field names. str1, str2, ... are filed names referring to property
      list name, field1, field2, ... are the corresponding values. Users are
      welcome to set any property as in the list in the constructor. If this
      is the case, users are also responsible for any confliction. In most
      cases, the following fields are required, :attr:`supercell`,
      :attr:`atoms`, :attr:`alist`, :attr:`xyzlist`, :attr:`ecut`. And the
      :attr:`name` field is recommended to be set. A typical example of the
      constructor can be found in :ref:`setup_molecule`.

   .. function:: Atom(asym)

      Constructor for Atom based on the atom symbol `asym`.

Class Property List
----------------------

.. class:: Molecule

   .. attribute:: name 

   .. attribute:: supercell

   .. attribute:: atoms

   .. attribute:: alist

   .. attribute:: xyzlist 

   .. attribute:: ecut

   .. attribute:: ecut2 

   .. attribute:: n1 

   .. attribute:: n2 

   .. attribute:: n3 

   .. attribute:: vext 

   .. attribute:: nspin 

   .. attribute:: temperature 

   .. attribute:: natoms 

   .. attribute:: vol 

   .. attribute:: nel 

