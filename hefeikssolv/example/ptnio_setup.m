%
% construct the SiH4 molecule (crystal)
%

%
% 1. construct atoms
%
a1 = Atom('Pt');
a2 = Atom('Ni');
a3 = Atom('O');
atomlist = zeros(1,9,'Atom');
for j = 1:6
   atomlist(j) = a1;
end
for j = 7:8
   atomlist(j) = a2;
end
atomlist(9) = a3;
%
% 2. set up primitive cell 
%
C = [
  0.19587405E+02   0.0000000E+00   0.0000000E+00
  0.0000000E+00   0.1066203E+02   0.0000000E+00
  0.0000000E+00   0.0000000E+00   0.9233592E+01
];
%
% 3. define the coordinates the atoms 
%
coefs = [
   0.648094   -0.017028029   -0.001622146
   0.646918    0.248561048    0.509308587
   0.652044    0.515916997   -0.007692607
   0.638433    0.755749261    0.503956325
   0.436396    0.241243102    0.829215792
   0.431923    0.506910875    0.330766318
   0.428427    0.000278279    0.332336185
   0.433327    0.750156053    0.840663815
   0.763344    0.248011715    0.163615671
];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','PtNiO', 'temperature', 3000 );
