%
%  construct bulk Silicon (crystal)
%
kssolvpptype('default')

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

autokpts = [ 6 6 6 0 0 0 ];

% 4. Configure the molecule (crystal)
cry = Crystal('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',20, 'name','Silicon', 'autokpts', autokpts );%'kpts',kpts, 'wks',wks );

opt = setksopt;
opt.verbose = 'off';
opt.maxscfiter = 99;

[cry,H,X,info] = scf(cry,opt);
