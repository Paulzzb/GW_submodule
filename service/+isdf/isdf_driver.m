% 
% License-Identifier: BSD-3-Clause
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/02/04
%
function isdf_driver(input_dir)
  % isdf_driver - driver for implementing ISDF\
  % Results are saved in ISDF_DB in input_dir
  % -----------------------------------------------------------------
  cleanup = output.push('isdf_driver');
  
  def = filename_map();

  % -----------------------------------------------------------------
  % Load GWinfo and config
  % 
  fName = fullfile(input_dir, def.GWinput);
  TEMP = load(fName);
  GWinfo = TEMP.GWgroundstate;
  fName = fullfile(input_dir, def.config);
  TEMP = load(fName);
  config = TEMP.config;
  TMPdir = config.CONTROL.storage_dir;
  clear TEMP;
  % -----------------------------------------------------------------
  if ~isempty(config.CONTROL.log_file)
    output.set_logfile(config.CONTROL.log_file);
  end
  % -----------------------------------------------------------------
  % Return if no isdf required
  if ~config.ISDF.isisdf
    output.msg('v0s', 'ISDF is not implemented!');
    output.msg('v0s', 'Return directly!');
    return
  end
  %
  % -----------------------------------------------------------------
  startisdf = tic;
  % Create descriptor
  dbroot = fullfile(TMPdir, def.isdf_database);
  new_describer = isdf_desc(GWinfo, config);
  % -----------------------------------------------------------------
  % Compare the descriptor, decide which isdf type to do
  isdf_flag_compute = [true, true, true, true]; % for type1, 2, 3
  db_dir = fullfile(TMPdir, def.isdf_database);
  if ~exist(db_dir, "dir") % No database for isdf, create one, then do the
                           % following calculation.
    output.msg('v0s', '%s', "No ISDF database found!\n");
    msg = sprintf("Create a new database at %s!\n", db_dir);
    output.msg('v0s', '%s', msg);
    meta = db_create(dbroot, new_describer, 1);
  else
    output.msg('v1s', '%s', "ISDF database found, check descriptor!\n");
    meta = db_read_meta(db_dir);
    old_describer = meta.desc;
    isdf_flag_compute = isdf_desc_compare(new_describer, old_describer);
  end
  meta.desc = new_describer;
  db_save(dbroot, meta);

  if isempty(GWinfo.psir) && config.ISDF.isisdf
    output.msg('v1s', '%s', 'Converting wavefunction from reciprocial space to real space ...');
    GWinfo.psir = get_wavefunc_real(GWinfo.psig, GWinfo.Ggrid4psig);
    output.msg('v2s', '%s', 'Wavefunction in real space prepared.');

  end
  %
  config.ISDF.dbroot = dbroot;
  dir = config.CONTROL.storage_dir;
  % GWgroundstate = GWinfo;
  % fName1 = fullfile(dir, def.GWinput);
  % save(fName1, 'GWgroundstate', '-v7.3', '-nocompression');
  fName2 = fullfile(dir, def.config);
  save(fName2, 'config');
  
  % Implement ISDF here!!!
  % Clause on is_helper, to decide if helper function output is needed.
  isdf_flag_compute(2) = false;
  indicater = [0, 0];
  starttimetype1 = tic;
  if isdf_flag_compute(1)
    [info1, info2, info3, info4] = isdf_sub("vc", indicater, dbroot, GWinfo, config); 
  end
  time1 = toc(starttimetype1);
  msg = sprintf('Time for computing vc: %f\n', time1);
  output.msg('v0s', '%s', msg);
  starttimetype2 = tic;
  if isdf_flag_compute(2)
    [info1, info2, info3, info4] = isdf_sub("vs", indicater, dbroot, GWinfo, config);
  end
  time2 = toc(starttimetype2);
  msg = sprintf('Time for computing vn: %f\n', time2);
  output.msg('v0s', '%s', msg);
  starttimetype3 = tic;
  if isdf_flag_compute(3)
    [info1, info2, info3, info4] = isdf_sub("ss", indicater, dbroot, GWinfo, config);
  end
  time3 = toc(starttimetype3);
  msg = sprintf('Time for computing nn: %f\n', time3);
  output.msg('v0s', '%s', msg);

  % Add some extra data

  if isdf_flag_compute(4)
    starttimetype4 = tic;
    hVh = isdf_vcVnn(dbroot,GWinfo,config);
    meta = db_read_meta(dbroot);
    meta = db_write(dbroot, meta, 4, "vcVnn", hVh);
    time4 = toc(starttimetype4);
    msg = sprintf('Time for computing vcVnn: %f\n', time4);
    output.msg('v0s', '%s', msg);
    msg = sprintf('ISDF driver done');
    output.msg('v0s', '%s', msg);
  end

  
  timeISDF = toc(startisdf);
  msg = sprintf('Time of ISDF driver: %f\n', timeISDF);
  output.msg('v0s', '%s', msg);
  %
end %function