function mol = finalize(mol)
% MOLECULE/FINALIZE Finalize function for molecule class
%    mol = FINALIZE(mol) returns a molecule class of with finalized fields.
%
%    See also Molecule.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if isempty(mol.name)
    mol.name = 'example';
end

if isempty(mol.supercell)
    error('supercell is must be set to start calculation !');
end

mol.vol = abs(det(mol.supercell));

if isempty(mol.atoms) || isempty(mol.alist)
    error('atomlist is must be set to start calculation !');
end

if size(mol.atoms,1) < size(mol.atoms,2)
    mol.atoms = mol.atoms(:);
end

if size(mol.alist,1) < size(mol.alist,2)
    mol.alist = mol.alist(:);
end

if ~isempty(mol.atoms)
    mol.natoms = sum(repmat(mol.alist,1,length(mol.atoms)) == ...
        repmat(unique(mol.alist)',length(mol.alist),1),1);
end

if isempty(mol.ecut)
    mol.ecut = 20.0;
end

if isempty(mol.ecut2)
    mol.ecut2 = 4*mol.ecut;
end

if isempty(mol.n1) || isempty(mol.n2) || isempty(mol.n3)
    C = mol.supercell;
    fackpt = sqrt(mol.ecut2*(2*meDef()))/pi;
    mol.n1 = ceil(fackpt*norm(C(1,:)));
    mol.n2 = ceil(fackpt*norm(C(2,:)));
    mol.n3 = ceil(fackpt*norm(C(3,:)));
end

if isempty(mol.vext)
    mol.vext = zeros(mol.n1,mol.n2,mol.n3);
end

if isempty(mol.smear)
    mol.smear = 'fd';
end

if isempty(mol.temperature)
    mol.temperature = 0.0;
end

if isempty(mol.alat)
    mol.alat = norm(mol.supercell(1,:));
end

% Initialize attributes related to electron spin
if isempty(mol.nspin)
    mol.nspin = 1;
end

if isempty(mol.tot_mag)
    mol.tot_mag = -1;
else
	error('tot_mag is not used by now and can not be set !');
end

if isempty(mol.lspinorb)
    mol.lspinorb = false;
end

if isempty(mol.noncolin)
    mol.noncolin = (mol.nspin == 4);
end

if mol.nspin == 4
    [mol.ux,mol.lsign] = compute_ux(mol.initmag);
end

if isempty(mol.domag)
    mol.domag = (mol.nspin == 2||mol.nspin == 4);
end

if isempty(mol.lsda)
    mol.lsda = (mol.nspin == 2);
end

if isempty(mol.ppvar)
    mol.ppvar = PpVariable(mol);
end

if isempty(mol.nel)
    mol.nel = sum(mol.natoms.*[mol.ppvar.venums]);
end

if ~isempty(mol.nbnd) && ~isempty(mol.extranbnd)
    error('nbnd and extranbnd can not be set together !');
end

if isempty(mol.funct)
    mol.funct = 'PZ';
end

if isempty(mol.xyzlist)
    error('xyzlist is must be set to start calculation !');
end

if isempty(mol.info)
    mol.info = '';
end

factor = setfunct(mol.funct);

if isempty(mol.omega)
    mol.omega = factor.omega;
end

if isempty(mol.exx_sr)
    mol.exx_sr = factor.exx_sr;
end

if isempty(mol.exx_lr)
    mol.exx_lr = factor.exx_lr;
end

if isempty(mol.dft_c)
    mol.dft_c = factor.dft_c;
end

end
