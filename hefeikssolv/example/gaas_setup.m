%
% construct the GaAs molecule (ethane)
%

%
% 1. construct atoms
%
a1 = Atom('Ga');
a2 = Atom('As');
atomlist = [a1 a1 a1 a1 a2 a2 a2 a2];
%
% 2. set up the supercell
%
BoxSize = 10.6826;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms 
%
xyzlist = [
0.0      0.0      0.0  
0.5      0.5      0.0  
0.5      0.0      0.5  
0.0      0.5      0.5  
0.25     0.25     0.25 
0.75     0.75     0.25 
0.75     0.25     0.75 
0.25     0.75     0.75 
]*C';

%
% 4. Configure the molecule
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','GaAs' );