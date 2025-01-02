function KSSOLV_Func_Tests()
clc;
clear;
clear class;
%rng(19890417);

tic;

run('../KSSOLV_startup.m');

kssolvpptype('default');

a1 = Atom('Si');
a2 = Atom('H');
alist(1) = a1;
alist(2) = a2;
alist(3) = a2;
alist(4) = a2;
alist(5) = a2;
%alist = [a1 a2 a2 a2 a2];
C = 10*eye(3);
coefs = [
 0.0     0.0      0.0
 0.161   0.161    0.161
-0.161  -0.161    0.161
 0.161  -0.161   -0.161
-0.161   0.161   -0.161
];
xyzmat = coefs*C';
mol = Molecule( 'supercell',C, 'atomlist',alist, 'xyzlist',xyzmat, ...
    'ecut',12.5, 'name','SiH4' );

[mol,H,X,~] = scf(mol);

pxyzlist = xyzmat + randn(size(xyzmat))/5;
molper = mol;
molper = set(molper,'xyzlist',pxyzlist);
relaxatoms(molper,1,H,X);
relaxatoms(molper,3,H,X);

dcm(mol);

trdcm(mol);

mol = Molecule( 'supercell',C, 'atomlist',alist, 'xyzlist',xyzmat, ...
    'ecut',12.5, 'name','SiH4', 'temperature', 1000 );

scf(mol);

mol2QEpw(mol);

autokpts = [ 2 2 2 0 0 0 ];

cry = Crystal( 'supercell',C, 'atomlist',alist, 'xyzlist',xyzmat, ...
    'ecut',12.5, 'name','SiH4k', 'autokpts', autokpts );

mol2QEpw(cry);

scf(cry);

cry = Crystal( 'supercell',C, 'atomlist',alist, 'xyzlist',xyzmat, ...
    'ecut',12.5, 'name','SiH4k', 'autokpts', autokpts, ...
    'temperature', 1000 );

mol2QEpw(cry);

scf(cry);

%%
kssolvpptype('pz-vbc');

a1 = Atom('Si');
a2 = Atom('H');
alist = [a1 a2 a2 a2 a2];
C = 10*eye(3);
coefs = [
 0.0     0.0      0.0
 0.161   0.161    0.161
-0.161  -0.161    0.161
 0.161  -0.161   -0.161
-0.161   0.161   -0.161
];
xyzmat = coefs*C';
mol = Molecule( 'supercell',C, 'atomlist',alist, 'xyzlist',xyzmat, ...
    'ecut',12.5, 'name','SiH4' );

scf(mol);

dcm(mol);

trdcm(mol);

mol = Molecule( 'supercell',C, 'atomlist',alist, 'xyzlist',xyzmat, ...
    'ecut',12.5, 'name','SiH4', 'temperature', 1000 );

scf(mol);

mol2QEpw(mol);

autokpts = [ 2 2 2 0 0 0 ];

cry = Crystal( 'supercell',C, 'atomlist',alist, 'xyzlist',xyzmat, ...
    'ecut',12.5, 'name','SiH4k', 'autokpts', autokpts );

mol2QEpw(cry);

scf(cry);

cry = Crystal( 'supercell',C, 'atomlist',alist, 'xyzlist',xyzmat, ...
    'ecut',12.5, 'name','SiH4k', 'autokpts', autokpts, ...
    'temperature', 1000 );

mol2QEpw(cry);

[cry,H,~,~] = scf(cry);

toc;

KSMolViewer(cry,H);

end
