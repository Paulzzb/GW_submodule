function KSSOLV_Tests()
clc;
clear;
clear class;
rng(19890417);

tic;

run('../KSSOLV_startup.m');

kssolvpptype('default');
tol = 1e-2;
checklistname = {'SiH4-SCF','SiH4-TRDCM', 'PtNiO-SCF','PtNiO-TRDCM', ...
    'Alanine-TRDCM','Al-SCF4M'};
checklist = zeros(length(checklistname),1);
opt = setksopt;
opt.verbose = 'off';
%%
% Check for SiH4 (Silane) molecule 
%
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

Eref = -6.1640047551846e+00;
[~,~,~,info] = scf(mol,opt);
checklist(1) = abs(info.Etotvec(end)-Eref)/abs(Eref)<tol;
check_report(checklistname{1},Eref,info.Etotvec(end),checklist(1));

Eref = -6.1639911648779e+00;
[~,~,~,info] = trdcm(mol,opt);
checklist(2) = abs(info.Etotvec(end)-Eref)/abs(Eref)<tol;
check_report(checklistname{2},Eref,info.Etotvec(end),checklist(2));

%%
% Check for PtNiO molecule 
%
a1 = Atom('Pt');
a2 = Atom('Ni');
a3 = Atom('O');
atomlist = [a1 a1 a1 a1 a1 a1 a2 a2 a3];
C = [
  0.19587405E+02   0.0000000E+00   0.0000000E+00
  0.0000000E+00   0.1066203E+02   0.0000000E+00
  0.0000000E+00   0.0000000E+00   0.9233592E+01
];
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
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','PtNiO' );

Eref = -1.9292720275602e+03;
[~,~,~,info] = scf(mol,opt);
checklist(3) = abs(info.Etotvec(end)-Eref)/abs(Eref)<tol;
check_report(checklistname{3},Eref,info.Etotvec(end),checklist(3));

Eref = -1.9459584306652e+03;
[~,~,~,info] = trdcm(mol,opt);
checklist(4) = abs(info.Etotvec(end)-Eref)/abs(Eref)<tol;
check_report(checklistname{4},Eref,info.Etotvec(end),checklist(4));

%%
% Check for Analine molecule 
%
a1 = Atom('H');
a2 = Atom('N');
a3 = Atom('C');
a4 = Atom('O');
atomlist = [a1, a2, a3, a3, a4, a3, a1, a4, a1, a1, a1, a1, a1];
C = [20 0 0; 0 15 0; 0 0 20];
xyzlist = [
  -0.133895    1.636840    1.067415
  -0.114950    0.674034    1.390795
  -0.129352   -0.247512    0.272425
   1.116503   -0.074057   -0.581443
   1.816148    0.910827   -0.622102
  -1.393741   -0.013479   -0.566277
   0.745649    0.588173    1.923262
   1.378429   -1.162243   -1.308285
   2.122853   -1.195908   -1.879225
  -0.127349   -1.285733    0.685615
  -1.453821   -0.733080   -1.415204
  -1.409682    1.017453   -0.989385
  -2.313398   -0.144089    0.049763
]/0.529177;
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist',xyzlist, ...
    'ecut',12.5, 'name','alanine');

Eref = -6.1124756054205e+01;
[~,~,~,info] = trdcm(mol,opt);
checklist(5) = abs(info.Etotvec(end)-Eref)/abs(Eref)<tol;
check_report(checklistname{5},Eref,info.Etotvec(end),checklist(5));

%%
% Check for Al molecule 
%
a = Atom('Al');
alist = [a a a a a a a a];
C = 20*eye(3);
xyzmat = [
0 0 0
10 0 0
0 10 0
0 0 10
10 10 10
10 10 0
0 10 10
10 0 10
];
mol = Molecule('supercell',C, 'atomlist',alist, 'xyzlist',xyzmat, ...
    'ecut',12.5, 'name','Aluminum cluster' );

mol = set(mol,'temperature',3000);
opts = setksopt;
opts.verbose = 'off';

Eref = -2.7593217153122e+02;
[~,~,~,info] = scf4m(mol,opts);
checklist(6) = abs(info.Etotvec(end)-Eref)/abs(Eref)<tol;
check_report(checklistname{6},Eref,info.Etotvec(end),checklist(6));

tottime = toc;

%%
fprintf('\n\n');
fprintf('=============================================================\n');
fprintf('                     Overall Report\n');
fprintf('-------------------------------------------------------------\n');
for it = 1:length(checklistname)
    fprintf('%18s   %d\n',checklistname{it},checklist(it));
end
fprintf('-------------------------------------------------------------\n');
fprintf('%18s   %f seconds\n','Total Running Time:',tottime);
fprintf('=============================================================\n');

%%
function check_report(testname,Eref,E,check)
fprintf('=============================================================\n');
fprintf('Test Name:        %s\n',testname);
fprintf('Reference Energy: %s\n',Eref);
fprintf('Current Energy:   %s\n',E);
fprintf('Energy Accuracy:  %e\n',abs(Eref-E)/abs(Eref));
if check
    fprintf('Test Result:      Passed\n');
else
    fprintf('Test Result:      Failed\n');
end
fprintf('=============================================================\n');
end

end
