function [sys,ebands,BH,BX] = eband(sys,options,kpts,rho,X)
% EBAND calculates the energy band.
%   PBE band: ebands = EBAND(sys,options,kpts,rho)
%   HSE band: ebands = EBAND(sys,options,kpts,rho,X)
%   calculates the smallest nbnd eigenvalues of a Hamiltonian
%   reconstructed by cry, rho and X at the given k-points.
%   ebands is a nkpts by nbnd matrix and each row of ebands is
%   the eigenvalues of the corresponding k-points.
%
%   rho and X should be provided explicitly,
%       the rho0 and X0 in options may be wrong.
% 
%   kpts format: {k_1x,k_1y,k_1z,num_kpts_k1_to_k2,k_1_name;
%                        k_2x,k_2y,k_2z,num_kpts_k2_to_k3,k_2_name}
%   kpts example: {0,0,0,21,'G';0.5,0,0,21,'X';0,0.5,0,21,'Y'}
%   nbnd: number of bands

if (nargin<5)
    X=[];
end

symkpts=reshape([kpts{:,1:3}],[],3);
nkpts=[kpts{:,4}];
kpoints=[kpts{1,1:3}];
for i = 1:size(kpts,1)-1
    nk=nkpts(i)-1;
    dk=(symkpts(i+1,1:3)-symkpts(i,1:3)).*((1:nk)'/nk);
    kpoints=[kpoints;symkpts(i,1:3)+dk];
end

%assert(~isempty(options.rho0),'Error: Non-SCF calculation must provide rho!');
options.rho0 = rho;
options.scfX = X;
options.X0 = [];
sys = set(sys,'kpts',kpoints,'scfkpts',sys.kpts);
sys = finalize(sys);

[sys,BH,BX,infokpts] = nscf4c(sys,options);

ebands = reshape(infokpts.Eigvals, sys.nbnd, []);

end
