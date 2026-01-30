% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/30
%
function isdf_driver(input_dir)
  % isdf_driver - driver for implementing ISDF\
  % Results are saved in ISDF_DB in input_dir
  % -----------------------------------------------------------------
  cleanup = QPlog_push('ISDF_driver');
  
  def = filename_map();
  % -----------------------------------------------------------------
  % Load GWinfo and config
  % 
  fName = fullfile(input_dir, def.GWinput);
  TEMP = load(fName);
  GWinfo = TEMP.GWgroundstate;
  config = TEMP.config;
  TMPdir = config.CONTROL.storage_dir;
  % -----------------------------------------------------------------
  % Return if no isdf required
  if ~config.ISDF.isisdf
    QPlog("ISDF is not implemented!\n", 0);
    QPlog("Return directly!\n", 0);
    return
  end
  %
  % -----------------------------------------------------------------
  % Create descriptor
  dbroot = fullfile(TMPdir, def.isdf_database);
  new_describer = isdf_desc(GWinfo, config);
  % -----------------------------------------------------------------
  % Compare the descriptor, decide which isdf type to do
  isdf_flag_compute = [true, true, true]; % for type1, 2, 3
  db_dir = fullfile(TMPdir, def.isdf_database);
  if ~exist(db_dir, "dir") % No database for isdf, create one, then do the
                           % following calculation.
    QPlog("No ISDF database found!\n", 0);
    msg = sprintf("Create a new database at %s!\n", db_dir);
    QPlog(msg, 0);
    meta = db_create(dbroot, new_describer, 1);
  else
    QPlog("ISDF database found, check descriptor!\n", 1);
    meta = db_read_meta(db_dir);
    old_describer = meta.desc;
    isdf_flag_compute = isdf_desc_compare(new_describer, old_describer);
  end

  % Implement ISDF here!!!
  QPlog('Converting wavefunction from reciprocial space to real space ...', 1);
  GWinfo.psir = get_wavefunc_real(GWinfo.psig, GWinfo.Ggrid4psig);
  QPlog('Wavefunction in real space prepared.', 2);
  
  % Clause on is_helper, to decide if helper function output is needed.
  indicater = [0, 0];
  if isdf_flag_compute(1)
    isdf_sub("vc", indicater, dbroot, GWinfo, config) 
  end
  if isdf_flag_compute(2)
    isdf_sub("vs", indicater, dbroot, GWinfo, config) 
  end
  if isdf_flag_compute(3)
    isdf_sub("ss", indicater, dbroot, GWinfo, config) 
  end
  %
  msg = sprintf('ISDF driver done');
  QPlog(msg, 0);
  %
end %function