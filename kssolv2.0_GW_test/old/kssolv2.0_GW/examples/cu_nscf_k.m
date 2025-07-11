%
% This script shows how to use KSSOLV to obtain the band structure of Cu
%

%
% first set up bulk Cu and compute its ground state
%
cu_scf_k;
%
% Then diagonalize the ground state Hamiltonian at several k-points
% and plot the band structure
%
endkpts = { 0.00  0.00  0.00  '\Gamma'
            0.00  0.50  0.50  'X'
            0.25  0.75  0.50  'W'
            0.50  0.50  0.50  'L'
            0.375 0.75  0.375 'K'};

plotband(cry,H.rho,info.efermi,15,endkpts,[0.4777 1]);
