% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/28
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
  end
  %
  % -----------------------------------------------------------------
  % Create descriptor
  new_describer = isdf_desc(GWinfo, config);
  % -----------------------------------------------------------------
  % Compare the descriptor, decide which isdf type to do
  isdf_flag_list = [false, false, false]; % for type1, 2, 3
  db_dir = fullfile(TMPdir, def.isdf_database);
  if ~exist(db_dir, "dir") % No database for isdf, create one, then do the
                           % following calculation.
    QPlog("No ISDF database found!\n", 0);
    msg = sprintf("Create a new database at %s!\n", db_dir);
    QPlog(msg, 0);
  else
    QPlog("ISDF database found, check descriptor!\n", 1);
    meta = db_read_meta(db_dir);
    old_describer = meta.describer;
    isdf_flag_list = isdf_desc_compare(new_describer, old_describer);
  end

  % Implement ISDF here!!!
  "Do isdf with main function, need to be done"
  quit

  % Save isdf result

end %function