%
% This script shows how to use the scf function to compute the ground state 
% of the bulk silicon using a small unit cell and k-point sampling
%
kssolvpptype('default')
%
%  construct the bulk Silicon (crystal) unit cell
%
%  1. construct atoms
a1 = Atom('Si');
atomlist = [a1 a1];

%  2. set up the primitive cell
C = 5.108*eye(3);

%  3. define the coordinates the atoms
xyzlist = [
    0.000000000         0.000000000         0.000000000
    0.600000000         0.500000000         0.500000000
    ]*C;

%
% Use 6x6x6 Monkhorst Pack k-point grid
%
autokpts = [ 6 6 6 0 0 0 ];

% 4. Configure the molecule (crystal) and set the kinetic energy cutoff
cry = Crystal('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecutwfc',20, 'name','Silicon', 'autokpts', autokpts );%'kpts',kpts, 'wks',wks );

%
% set the options by default and modify some of the options to be passed to the scf function
opt = setksopt;
opt.verbose = 'off';
opt.maxscfiter = 99;
%
% call scf to compute the ground state
[cry,H,X,info] = scf(cry,opt);
