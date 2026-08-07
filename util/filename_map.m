% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/08/04 ZZ

function def = filename_map()
  def.kssolvinput = 'groundstate.mat';
  def.qeinput = 'still_undergoing';
  def.data = 'data.mat';
  def.GWinput = 'GWinput.mat';
  def.config = 'config.mat';
  def.stage = 'relay_stage.mat';
  def.isdfvc = 'isdf_typevc.mat';
  def.isdfvs = 'isdf_typevs.mat';
  def.isdfss = 'isdf_typess.mat';
  def.QPenergy = 'QPenergy';
  def.isdf_database = 'ISDF_DB';
  def.isdf_report_dir = 'isdf_report';
  def.cond_report = 'o-ISDF_cond';
  def.hf_report = 'o-ISDF_HF_id%d';
  def.adaptive_report = 'o-ISDF_adaptive_id%d';
end % EOF
