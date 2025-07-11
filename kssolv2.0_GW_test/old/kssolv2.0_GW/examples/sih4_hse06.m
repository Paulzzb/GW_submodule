%
% This example shows how to call scf to perform a hybrid functional KSDFT calculation for 
% the SiH4 molecule using the HSE06 hybrid functional.
%

% set up the SiH4 molecule and 
% define some options to be used in scf
%
sih4_setup;
%
% specify the hybrid functional to be used for KSDFT calculation
mol.funct = 'HSE06';
%
% define the corresponding pseudopotential to be used
ppvar = PpVariable(mol);
ppvar.funct = 'HSE06';
mol.ppvar = ppvar;

%
% modify some of the options to be passed to the scf function 
% 
options = setksopt('betamix',0.6,'mixdim',8,'maxeigsiter',10,'scftol',1e-5,...
                    'eigstol',1e-6,'maxscfiter',20,'maxcgiter',10,'maxphiiter',10);

% specify using the Adaptive Compressive Exchange operator to 
% reduce the cost of the SCF calculation.
options.useace = 1;

%
% call the scf function to perform the SCF calcuation
[mol,H,X,info] = scf(mol, options);
