function sys = kssolv2GW_initsys(mol)
  % kssol
  F = KSFFT(mol);
  [sys.ng, sys.nr] = size(F);
  clear F
  sys.vol = mol.vol;
  sys.xyzlist = mol.xyzlist;
  sys.n1 = mol.n1;
  sys.n2 = mol.n2;
  sys.n3 = mol.n3;
  sys.ne = mol.nel;
  sys.supercell = mol.supercell;

end % EOF
