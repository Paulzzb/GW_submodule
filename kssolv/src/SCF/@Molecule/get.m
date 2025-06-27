function val = get(mol,attr_name)
% MOLECULE/GET Obtain attribute of molecule class
%    val = GET(mol,attr_name) returns attribute in mol
%
%    See also Molecule.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

switch attr_name
  case 'name'
    val = mol.name;
  case 'supercell'
    val = mol.supercell;
  case 'vol'
    val = mol.vol;
  case 'atoms'
    val = mol.atoms;
  case 'alist'
    val = mol.alist;
  case 'natoms'
    val = mol.natoms;
  case 'ecut'
    val = mol.ecut;
  case 'ecut2'  
    val = mol.ecut2;
  case 'n1'
    val = mol.n1; 
  case 'n2'
    val = mol.n2;
  case 'n3'
    val = mol.n3;
  case 'vext'
    val = mol.vext;
  case 'smear'
    val = mol.smear;
  case 'temperature'
    val = mol.temperature;
  case 'alat'
    val = mol.alat;
  case 'nspin'
    val = mol.nspin;
  case 'lspinorb'
    val = mol.lspinorb;
  otherwise 
    error('invalid attribute requested for Molecule');
end;

