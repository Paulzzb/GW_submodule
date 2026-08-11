% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/05/06 ZZ

function data = load_qe_groundstate(dirin)
% load_qe_groundstate - Load groundstate data from a QE pw.x output.
% User must have previously run a QE pw.x calculation.
myneed = load_qe_from_folder(dirin);
% error('load_qe_groundstate not implemented yet.');
data = struct();
data.rhor = myneed.rhor;
data.Vxc = myneed.vxc;
data.ev = myneed.ev;
data.psig = myneed.psig;
data.sys = myneed.sys;
data.occupation = myneed.occupation;
reciprocal_grid_info = struct();
reciprocal_grid_info.fftgrid = [myneed.n1, myneed.n2, myneed.n3];
reciprocal_grid_info.vol = myneed.vol;
reciprocal_grid_info.idxnz = myneed.idxnz;
reciprocal_grid_info.wfncut = myneed.wfncut;
reciprocal_grid_info.xyz = myneed.mill;

data.reciprocal_grid_info = reciprocal_grid_info;
data.nkibz = myneed.nkibz;
data.kibz = myneed.kibz;
data.kweight = myneed.kweight;
data.nspin = myneed.nspin;
data.nspinor = myneed.nspinor;

symm = myneed.syms;
data.syms = symm;

if isfield(myneed, 'xyz')
  data.xyz = myneed.xyz;
end
if isfield(myneed, 'atom_symbol')
  data.atom_symbol = myneed.atom_symbol;
end

end

