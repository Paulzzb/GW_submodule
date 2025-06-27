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
cell_au = 10; % cell_A is too small, only for test
C = cell_au*eye(3);
%
% 3. Define the coordinates of atoms 
%
bond = 1.5/cell_au;
coefs = [
   0.5000000000000000    0.5000000000000000    0.5-bond/2
   0.5000000000000000    0.5000000000000000    0.5+bond/2
];
xyzlist = coefs*C';

funct = {'PZ','PBE','HSE06'};
method = {1,3,4,5};
refer = {
    {1.46376357, 1.46389002, 1.46376357, 1.46394053},...
    {1.43403984, 1.43402226, 1.43403984, 1.43399044},...
    {1.42560203, 1.42570375, 1.42575101, 1.42565160}};
check = zeros(size(funct,2)*size(method,2),5);
for fi = 1:size(funct,2)
    for ni = 1:size(method,2)
        %
        % 4. Configure the molecule (crystal)
        %
        mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist',xyzlist, ...
            'ecut',20, 'name','H2', 'funct',funct{fi}, 'n1',45, 'n2',45, 'n3',45, ...
            'extranbnd',0, 'nspin',1);
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
        
        opt.relaxmethod = method{ni};
        opt.relaxtol = 1e-4;
        opt.maxrelaxiter = 200;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % measure the initial bond length between two H atoms
        bondlength = norm(mol.xyzlist(1,:)-mol.xyzlist(2,:));
        fprintf('starting bond length = %11.8e (Bohr)\n', bondlength);
        
        % perform geometry optimization of the bond length between two H atoms
        [molopt,Hopt,Xopt,infoopt] = relaxatoms(mol,opt);
        relaxHistory = infoopt.relaxHistory; % The information of each optimization step can be found in infoopt.relaxHistory
        
        % measure the optimized bond length between two H atoms
        bondlength = norm(molopt.xyzlist(1,:)-molopt.xyzlist(2,:));
        fprintf('optimized bond length = %11.8e (Bohr)\n', bondlength);
        
        check((fi-1)*size(method,2)+ni,1) = bondlength;
        check((fi-1)*size(method,2)+ni,2) = refer{fi}{ni};
        check((fi-1)*size(method,2)+ni,3) = check((fi-1)*size(method,2)+ni,1) - check((fi-1)*size(method,2)+ni,2);
        check((fi-1)*size(method,2)+ni,4) = abs(check((fi-1)*size(method,2)+ni,3)) < 1E-8;
        check((fi-1)*size(method,2)+ni,5) = infoopt.converge;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 无初始化
% PZ    1E-4    KSSOLV: 1.46376357(1)   1.46389002(3)   1.46376357(4)   1.46394053(5-0.05fs)    1.46357537(5-0.1fs)   QE: 1.46375146 (bfgs)
%       1E-10           1.46374688(1)*  1.46374612(3)*  1.46374692(4)*  1.46374677(5-0.05fs)    1.46374677(5-0.1fs)   
% PBE   1E-4    KSSOLV: 1.43403984(1)   1.43402226(3)   1.43403984(4)   1.43399044(5-0.05fs)    1.43390109(5-0.1fs)   QE: 1.43397663 (bfgs)
% HSE06 1E-4    KSSOLV: 1.42560203(1)   1.42570375(3)   1.42575101(4)   1.42565160(5-0.05fs)    1.42557463(5-0.1fs)   QE: 1.42568541 (bfgs)
%       1E-6            1.42560203(1)   1.42560164(3)   1.42560203(4)   1.42560283(5-0.05fs)    1.42560127(5-0.1fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 有初始化（原版）——初猜算两遍
% PZ    1E-4    KSSOLV: 1.46376366(1)   1.46389019(3)   1.46376366(4)   QE: 1.46375146 (bfgs)
% PBE   1E-4    KSSOLV: 1.43404010(1)   1.43402239(3)   1.43404010(4)   QE: 1.43397663 (bfgs)
% HSE06 1E-4    KSSOLV: 1.42575170(1)*  1.42734747(3)*  1.42575170(4)   QE: 1.42568541 (bfgs)   *表示无法收敛
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%