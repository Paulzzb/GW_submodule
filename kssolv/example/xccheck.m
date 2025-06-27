clc;
clear all;
close all;
%
% 0. Choose pseudopotential
%
kssolvpptype('default');
%
% 1. Construct atoms
%
a1 = Atom('H');
atomlist = [a1, a1];
%
% 2. Set up supercell
%
cell_au = 10;
C = cell_au*eye(3);
%
% 3. Define the coordinates of atoms 
%
bond = 1.4011/cell_au;
coefs = [
   0.5000000000000000    0.5000000000000000    0.5-bond/2
   0.5000000000000000    0.5000000000000000    0.5+bond/2
];
xyzlist = coefs*C';

funct = {'PZ','PW','PBE','HSE06','PBE0','HF','LC-WPBE','LC-PBE0'};
nspin = {1,2,4};
refer = {
    {-2.2614644421165e+00, -2.2614644421695e+00, -2.2614644421922e+00},...
    {-2.2607795623778e+00, -2.2607795624214e+00, -2.2607795624556e+00},...
    {-2.3180802593943e+00, -2.3180845190847e+00, -2.3180845189553e+00},...
    {-2.3268843157874e+00, -2.3255975632283e+00, -2.3255975690828e+00},...
    {-2.3293948627558e+00, -2.3293956137584e+00, -2.3293956157233e+00},...
    {-2.2761103313358e+00, -2.2761003709775e+00, -2.2761003805482e+00},...
    {-2.3542114105516e+00, -2.3576983060113e+00, -2.3576983152834e+00},...  % exxgkk的长程部分已产生变化（差距-0.0084）
    {-2.3494384919709e+00, -2.3524105402346e+00, -2.3524105406662e+00}};    % exxgkk的长程部分已产生变化（差距-0.0084）
check = zeros(size(funct,2)*size(nspin,2),4);

for fi = 1:size(funct,2)
    for ni = 1:size(nspin,2)
    %
    % 4. Configure the molecule (crystal)
    %
    mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist',xyzlist, ...
        'ecut',15, 'name','H2', 'funct',funct{fi}, 'n1',36, 'n2',36, 'n3',36, ...
        'extranbnd',0, 'nspin',nspin{ni});
    %
    % 5. Configure the options
    %
    opt = setksopt();

    opt.eigmethod = 'davidson_qe';
    opt.scftol    = 1e-9;
    opt.phitol    = 1e-9;
    opt.maxscfiter= 300;
    opt.maxphiiter= 300;
    opt.mixtype   = 'broyden_modify';
    opt.what2mix  = 'rho';
    opt.betamix   = 0.7;
    opt.mixdim    = 8;
    opt.brank     = 8;

    opt.exxmethod = 'ace';
    opt.x_gamma_extrapolation = false;

    [mol,H0,X0,info] = scf(mol,opt);
    check((fi-1)*size(nspin,2)+ni,1) = info.Etot*2;
    check((fi-1)*size(nspin,2)+ni,2) = refer{fi}{ni};
    check((fi-1)*size(nspin,2)+ni,3) = check((fi-1)*size(nspin,2)+ni,1) - check((fi-1)*size(nspin,2)+ni,2);
    check((fi-1)*size(nspin,2)+ni,4) = abs(check((fi-1)*size(nspin,2)+ni,3)) < 1E-10; % 1E-12
    end
end
%PZ:        -2.2614644421165e+00@abs-2.2614644421695e+00(nspin=2)   -2.2614644421922e+00(nspin=4)
%           -2.26146444             -2.26146444
%PW:        -2.2607795623770e+00@abs-2.2607795624214e+00(nspin=2)   -2.2607795624556e+00(nspin=4) @self
%           -2.26077956             -2.26077956
%PBE:       -2.3180802593943e+00    -2.3180845190847e+00(nspin=2)   -2.3180845189553e+00(nspin=4)
%           -2.31808025             -2.31808452
%HSE06:     -2.3268843157874e+00    -2.3255975632283e+00(nspin=2)   -2.3255975690828e+00(nspin=4)
%           -2.32688431             -2.32559336
%PBE0:      -2.3293948627558e+00    -2.3293956137584e+00(nspin=2)   -2.3293956157233e+00(nspin=4)
%           -2.32939486             -2.32939446
%HF:        -2.2761103313358e+00    -2.2761003709775e+00(nspin=2)   -2.2761003805482e+00(nspin=4)
%           -2.27611033             -2.27611033
%LC-WPBE:   -2.3542114105516e+00    -2.3576983060113e+00(nspin=2)   -2.3576983152834e+00(nspin=4)
%LC-PBE0:   -2.3494384919709e+00    -2.3524105402346e+00(nspin=2)   -2.3524105406662e+00(nspin=4)