%
%  construct bulk Silicon (crystal)
%
kssolvpptype('default')

%  1. construct atoms
a1 = Atom('Si');
atomlist = [a1 a1 a1 a1 a1 a1 a1 a1];

%  2. set up the primitive cell
C = 10.216*eye(3);

%  3. define the coordinates the atoms
xyzlist = [
    0.000000000         0.000000000         0.000000000
    0.000000000         0.500000000         0.500000000
    0.500000000         0.000000000         0.500000000
    0.500000000         0.500000000         0.000000000
    0.750000000         0.250000000         0.750000000
    0.250000000         0.250000000         0.250000000
    0.250000000         0.750000000         0.750000000
    0.750000000         0.750000000         0.250000000
    ]*C;

autokpts = [ 2 2 2 0 0 0 ];

% 4. Configure the molecule (crystal)
cry = Crystal('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',30, 'name','Silicon', 'autokpts', autokpts );%'kpts',kpts, 'wks',wks );

opt = setksopt;
opt.verbose = 'off';
opt.maxscfiter = 10;

% QE ref: Etot = -63.48018027
[cry,H,X,info] = scf(cry,opt);
