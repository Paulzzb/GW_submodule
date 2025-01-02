% 1. Set up atoms for the molecule
a1 = Atom('H');
a2 = Atom('Si');
atomlist = [a2 a1 a1 a1 a1];

% 2. Set up the supercell in the unit of Bhor
C = 10*eye(3);

% 3. Set up the energy cut for frequency domain in the unit of Hartree
ecut = 12.5;

% 4. Set the coordinates for each atom 
xyzlist = [  0.00   0.00    0.00
             1.61   1.61    1.61
            -1.61  -1.61    1.61
             1.61  -1.61   -1.61
            -1.61   1.61   -1.61 ];

% 5. Set the molecule SiH4 from values defined above
mol = Molecule( 'name','SiH4', 'supercell',C, 'atomlist',atomlist, ...
    'xyzlist' ,xyzlist, 'ecut',ecut );
