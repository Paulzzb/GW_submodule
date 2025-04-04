function [mol,H,X,Hprec,varargout] = iterinit(mol,varargin)
% ITERINIT generates initial values for the iterative method.
%    [mol,H,X,Hprec] = ITERINIT(mol) generates the initial Hamiltonian,
%    wave functions and preconditioner of the corresponding Hamiltonian for
%    molecule mol.
%
%    [mol,H,X,Hprec] = ITERINIT(mol,rho) generates the initial Hamiltonian,
%    wave functions and preconditioner of the corresponding Hamiltonian for
%    molecule mol together with density function rho.
%
%    [mol,H,X,Hprec] = ITERINIT(mol,X) generates the initial Hamiltonian,
%    wave functions and preconditioner of the corresponding Hamiltonian for
%    molecule mol together with wave function X.
%
%    [mol,H,X,Hprec] = ITERINIT(mol,rho,X) generates the initial
%    Hamiltonian, wave functions and preconditioner of the corresponding
%    Hamiltonian for molecule mol together with density function rho and
%    wave function X.
%
%    [mol,H,X,Hprec] = ITERINIT(mol,rho,X,ncol) generates X with ncol
%    columns.
%
%    [cry,BH,BX,BHprec] = ITERINIT(cry) generates the initial Bloch
%    Hamiltonian, Bloch wave functions and the corresponding preconditioner
%    for the given crystal cry.
%
%    [cry,BH,BX,BHprec] = ITERINIT(cry,rho) generates the initial Bloch
%    Hamiltonian, Bloch wave functions and the corresponding preconditioner
%    for the given crystal cry together with density function rho.
%
%    [cry,BH,BX,BHprec] = ITERINIT(cry,BX) generates the initial Bloch
%    Hamiltonian, Bloch wave functions and the corresponding preconditioner
%    for the given crystal cry together with Bloch wave function BX.
%
%    [cry,BH,BX,BHprec] = ITERINIT(cry,rho,BX) generates the initial Bloch
%    Hamiltonian, Bloch wave functions and the corresponding preconditioner
%    for the given crystal cry together with density function rho and Bloch
%    wave function BX.
%
%    [cry,BH,BX,BHprec] = ITERINIT(cry,rho,BX,ncols) generates BX with each
%    number in ncols being the number of columns for each k-points.
%
%    [...,nocc] = ITERINIT(...) returns the total number of occupied
%    colmuns for wave function or Bloch wave function together with other
%    inputs and outputs.
%
%   See also Ham, Wavefun, BlochHam, BlochWavefun, genX0, genprec.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if isempty(mol.ppvar)
    mol.ppvar  = PpVariable(mol);
    ne = 0;
     for it = 1:length(mol.atoms)
        ne = ne + mol.ppvar.venums(it)*mol.natoms(it);
    end
    if ne ~= mol.nel
        warning(['The number of valence electrons in Pp file is ' ...
            'different from Periodic Table']);
        mol = set(mol,'nel',ne);
    end
end

% Set the number of up and down electrons,only useful in constrained 
% magnetization optimization which has not been realized by now
mol = set_nelup_neldw(mol);

if mol.temperature > 0
    tmul = 1;
else
    tmul = 0;
end

%number of valence bands
if isa(mol,'Crystal')
    if ~mol.noncolin
        nocc = mol.nel/2*ones(mol.nkpts*mol.nspin,1);
        ncol = max(mol.nbnd, ceil(nocc + tmul*round(max(4,0.2*mol.nel/2))));
    else
        nocc = mol.nel*ones(mol.nkpts,1);
        ncol = max(mol.nbnd, ceil(nocc + tmul*round(max(4,0.2*mol.nel))));
    end

    if nargin == 1
        rho = [];
        X = [];
    elseif nargin == 2
        if isa(varargin{1},'BlochWavefun')
            rho = [];
            X = varargin{1};
        else
            rho = varargin{1};
            X = [];
        end
    elseif nargin == 3
        rho = varargin{1};
        X = varargin{2};
    else
        rho = varargin{1};
        X = varargin{2};
        ncol = ceil(varargin{3});
    end
    
    if isempty(rho)
        H = BlochHam(mol);
    else
        H = BlochHam(mol,rho);
    end
else
    if ~mol.noncolin
        nocc = mol.nel/2*mol.nspin;
        ncol = max(mol.nbnd, ceil(nocc/mol.nspin + tmul*round(max(4,0.2*mol.nel/2))));
    else
        nocc = mol.nel;
        ncol = max(mol.nbnd, ceil(nocc + tmul*round(max(4,0.2*mol.nel))));
    end

    if nargin == 1
        rho = [];
        X = [];
    elseif nargin == 2
        if isa(varargin{1},'Wavefun')
            rho = [];
            X = varargin{1};
        else
            rho = varargin{1};
            X = [];
        end
    elseif nargin == 3
        rho = varargin{1};
        X = varargin{2};
    else
        rho = varargin{1};
        X = varargin{2};
        ncol = ceil(varargin{3});
    end
    
    if isempty(rho)
        H = Ham(mol);
    else
        H = Ham(mol,rho);
    end
end

if ~mol.abnd
    mol = set(mol,'nbnd',max(ncol));
end

if isempty(X) 
    %X = genX0(mol,ncol);
    [X, H] = genX0_random(mol,H);
else
    if ~mol.lsda || isa(X, 'BlochWavefun')
        if any(ncols(X) < nocc + tmul*round(max(4,0.2*mol.nel/2*(1+mol.noncolin))))
            error('The given initial X is not valid.');
        end
    else
        if any(ncols(X{1}) < nocc/2 + tmul*round(max(4,0.2*mol.nel/2)))...
           || any(ncols(X{2}) < nocc/2 + tmul*round(max(4,0.2*mol.nel/2)))
            error('The given initial X is not valid.');
        end
    end 
end

if ~mol.noncolin
    Hprec = genprec(H);
else
    Hprec = genprec(H);
    if isa( H, 'BlochHam' )
        for ik = 1:H.nkpts
            Hprec{ik}.psi = repmat(Hprec{ik}.psi,2,1);
        end
    else
        Hprec.psi = repmat(Hprec.psi,2,1);
    end
end

if nargout == 5
    varargout{1} = sumel(nocc);
end

end
