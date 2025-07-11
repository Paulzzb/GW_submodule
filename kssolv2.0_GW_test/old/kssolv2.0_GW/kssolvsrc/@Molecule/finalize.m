function mol = finalize(mol)
% MOLECULE/FINALIZE Finalize function for molecule class
%    mol = FINALIZE(mol) returns a molecule class of with finalized fields.
%
%    See also Molecule.

if isempty(mol.ecutrho)
    mol.ecutrho = 4*mol.ecutwfc;
end

if isempty(mol.n1) || isempty(mol.n2) || isempty(mol.n3)
    C = mol.supercell;
    fackpt = sqrt(mol.ecutrho*(2*meDef()))/pi;
    mol.n1 = ceil(fackpt*norm(C(1,:)));
    mol.n2 = ceil(fackpt*norm(C(2,:)));
    mol.n3 = ceil(fackpt*norm(C(3,:)));
end

mol.gridwfc = Ggrid(mol, mol.ecutwfc);
mol.gridrho = Ggrid(mol, mol.ecutrho);

if size(mol.atoms,1) < size(mol.atoms,2)
    mol.atoms = mol.atoms(:);
end

if size(mol.alist,1) < size(mol.alist,2)
    mol.alist = mol.alist(:);
end

if isempty(mol.vext)
    mol.vext = zeros(mol.n1,mol.n2,mol.n3);
end

if isempty(mol.nspin)
    mol.nspin = 1;
end

if isempty(mol.temperature)
    mol.temperature = 0.0;
end

if ~isempty(mol.atoms)
    mol.natoms = sum(repmat(mol.alist,1,length(mol.atoms)) == ...
        repmat(unique(mol.alist)',length(mol.alist),1),1);
end

mol.vol = abs(det(mol.supercell));

if isempty(mol.nel) && ~isempty(mol.atoms)
    mol.nel = sum(mol.natoms.*[mol.atoms.venum]);
end

end
