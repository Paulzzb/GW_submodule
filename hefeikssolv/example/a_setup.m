%
% construct the SiH4 (Silane) molecule 
%
%kssolvpptype('default');
kssolvpptype('pz-hgh', 'UPF');
%
% 1. construct atoms
%
a1 = Atom('H');
a2 = Atom('C');
a3 = Atom('N');
atomlist = [ ...
a1, a1, a1, a1, a1, ...
a2, a2, a2, a2, a2, ...
a3, a3, a3, a3, a3 ...
];
%
% 2. set up supercell
%
C = 28.345891*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
     0.418949991         0.301939994         0.506910026
     0.702019989         0.397709996         0.502229989
     0.697579980         0.515619993         0.498519987
     0.456539989         0.728969991         0.492049992
     0.330819994         0.613150001         0.496650010
     0.448139995         0.368990004         0.504610002
     0.432610005         0.514540017         0.499920011
     0.524469972         0.528439999         0.499159992
     0.577440023         0.451469988         0.501380026
     0.466949999         0.657060027         0.494619995
     0.390859991         0.436560005         0.502640009
     0.536970019         0.372269988         0.504119992
     0.666649997         0.455419987         0.500729978
     0.544629991         0.617340028         0.495860010
     0.397069991         0.598460019         0.497000009
];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',10.0, 'name','DNA' );