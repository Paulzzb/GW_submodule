cd ../../../
KSSOLV_startup
cd silicon
si2
cd ../src/GW
gvec0 = gvec();
gvec1 = gvec(mol);
% gvec2 = gvec(mol, 'method', 'bgw');
gvec3 = gvec(mol, 'ecut', 14);

