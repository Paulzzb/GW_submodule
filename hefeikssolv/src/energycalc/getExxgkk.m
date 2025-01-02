function options = getExxgkk(mol,options)
% Calculate Exxgkk with k points, it calls getExxgkk_q.
options.exxcut=mol.ecut2;
grid = Ggrid(mol, options.exxcut);
idxnz=grid.idxnz;
n1 = mol.n1;n2 = mol.n2;n3 = mol.n3;
n123 = n1*n2*n3;

if isa(mol,'Crystal')
    if (~options.aceconv && ~options.store_tensors) || ...
           strcmp(options.exxmethod,'normal') || ...
               strcmp(options.exxmethod,'ace') 
        nk_psi = mol.nkpts;
        kpts_psi = mol.kpts;
        if isfield(options,'scfX') && ~isempty(options.scfX)
            nk_phi = options.scfX.nkpts;
            kpts_phi = mol.scfkpts;
        else
            nk_phi = nk_psi;
            kpts_phi = kpts_psi;
        end
        exxgkk=zeros(n123,nk_psi,nk_phi);
        for i=1:nk_psi
            for j=1:nk_phi
                %In QE, k+q=kq, i is index of k, j is index of kq.
                %The input argument of the formulation of exxgkk is k-kq=-q.
                %See function g2_convolution of exx_base.f90 in QE.
                q=kpts_psi(i,:)-kpts_phi(j,:);
                [exxgkk(idxnz,i,j),options]=getExxgkk_q(mol,options,q);
            end
        end
    else
        [kpts2,exxidxnz]=kgrid(mol);
        nkpts2=size(kpts2,1);
        exxgkk=zeros(n123,nkpts2);
        for ik=1:nkpts2
            [exxgkk(idxnz,ik),options]=getExxgkk_q(mol,options,kpts2(ik,:));
        end
        options.exxidxnz=exxidxnz;
    end
else
    exxgkk = zeros(n123,1);
    [exxgkk(idxnz), options] = getExxgkk_q(mol,options);
end
options.exxgkk=exxgkk;
end

function [exxgkk, options] = getExxgkk_q(mol,options,q)
%    Compute eigenvalues of the exact exchange operator
%    Equivalent to the function "g2_convolution" in exx_base.f90 of QE.
% Input:
%    mol  --- Molecule information
%
% Output:
%    Exxgkk --- Eigenvalues of the exact change operator for solving Poisson
%               equations

% Gygi-Baldereschi regularization
Mu = 0.106;     % Screening Parameter
epsDiv = 1e-8;  % Divergence threshold
grid = Ggrid(mol, options.exxcut);

if isprop(mol,'nqs') && ~isempty(mol.nqs)
    nq=mol.nqs;
elseif isa(mol,'Crystal')
    nq=mol.nkxyz;
else
    nq=[1,1,1];
end
if options.x_gamma_extrapolation
    grid_factor=8/7;
else
    grid_factor=1;
end
[exxDiv,options]=getexxDiv(mol,options);
if (nargin < 3)
    gkxyz=[grid.gkx,grid.gky,grid.gkz];
    gkk = grid.gkk;
else
    gkxyz=[grid.gkx,grid.gky,grid.gkz]+q;
    gkk = sum(gkxyz.^2,2);
end
grid_factor_track = ones(size(gkk));
if options.x_gamma_extrapolation
    x = gkxyz*(mol.supercell.'.*nq/pi/4);
    on_double_grid=sum(abs(x-round(x)),2)<3*eps;
    grid_factor_track(on_double_grid)=0;
    grid_factor_track(~on_double_grid)=grid_factor;
end
idx = gkk < epsDiv;
exxgkk = 4 * pi ./ gkk .* (1 - exp(-gkk / (4*Mu^2))).*grid_factor_track;
exxgkk(idx) = -exxDiv;
if ~options.x_gamma_extrapolation
    exxgkk(idx) = exxgkk(idx) + pi / Mu^2;
end

    function [exxDiv,options] = getexxDiv(mol,options)
        %Calculate the divergent G=0 term.
        %Equivalent to the function "exx_divergence" in exx_base.f90 of QE.
        if isfield(options,'exxDiv')
            exxDiv=options.exxDiv;
        else
            % Gygi-Baldereschi regularization
            exxAlpha = 10 / options.exxcut*2;
            rcell=inv(mol.supercell)'*2*pi;
            exxDiv=0;
            for iq1=0:nq(1)-1
                for iq2=0:nq(2)-1
                    for iq3=0:nq(3)-1
                        qxyz=[iq1/nq(1),iq2/nq(2),iq3/nq(3)]*rcell;
                        gkxyz=[grid.gkx,grid.gky,grid.gkz]+qxyz;
                        if options.x_gamma_extrapolation
                            x = gkxyz*(mol.supercell.'.*nq/pi/4);
                            on_double_grid=sum(abs(x-round(x)),2)<3*eps;
                        else
                            on_double_grid = zeros(size(gkxyz,1),1);
                        end
                        gkk = sum(gkxyz.*gkxyz,2);
                        gkki = exp(-exxAlpha * gkk) ./ gkk .* (1 - exp(-gkk / (4*Mu^2)));
                        idx = (gkk > epsDiv)&~on_double_grid;
                        exxDiv=exxDiv+sum(gkki(idx))*grid_factor;
                    end
                end
            end
            if ~options.x_gamma_extrapolation
                exxDiv = exxDiv + 1/(4*Mu^2);
            end
            exxDiv = 4 * pi * exxDiv;
            nqq  = 1e5; % What is this?
            dq = 5 / sqrt(exxAlpha) / nqq;
            qt2 = (((1:nqq+1)' - 0.5) * dq).^2;
            aa = -sum(exp(-(exxAlpha + 1/(4*Mu^2)) * qt2)) * dq * 2 / pi + 1/sqrt(exxAlpha * pi);
            exxDiv = exxDiv - mol.vol * aa*prod(nq);
            options.exxDiv=exxDiv;
        end
    end
end
