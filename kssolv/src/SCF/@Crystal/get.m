function val = get(cry,attr_name)
% CRYSTAL/GET Obtain attribute of Crystal class
% 	val = GET(cry,attr_name) returns attribute in cry
%
%   See also Crystal.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

switch attr_name
  case 'name'
    val = cry.name;
  case 'supercell'
    val = cry.supercell;
  case 'vol'
    val = cry.vol;
  case 'atoms'
    val = cry.atoms;
  case 'alist'
    val = cry.alist;
  case 'natoms'
    val = cry.natoms;
  case 'ecut'
    val = cry.ecut;
  case 'ecut2'
    val = cry.ecut2;
  case 'n1'
    val = cry.n1;
  case 'n2'
    val = cry.n2;
  case 'n3'
    val = cry.n3;
  case 'vext'
    val = cry.vext;
  case 'smear'
    val = cry.smear;
  case 'temperature'
    val = cry.temperature;
  case 'alat'
    val = cry.alat;
  case 'nspin'
    val = cry.nspin;
  case 'lspinorb'
    val = cry.lspinorb;
  case 'kpts'
    val = cry.kpts;
  case 'nkpts'
    val = cry.nkpts;
  case 'wks'
    val = cry.wks;
  case 'nkxyz'
    val = cry.nkxyz;
  case 'scfkpts'
    val = cry.scfkpts;
  case 'nqs'
    val = cry.nqs;
  otherwise 
    error('invalid attribute requested for Crystal');
end
