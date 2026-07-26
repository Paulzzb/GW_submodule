% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/23

function info = save2db(filePath, stage)
% Save staged relay snapshot to file in MAT format (-v7.3 nocompression).
% If stage is omitted, use the staged cache from appdata.

  if nargin < 1 || isempty(filePath)
    error('relay::save2db requires filePath.');
  end

  if nargin < 2 || isempty(stage)
    if ~isappdata(0, 'GW_RELAY_STAGE')
      error('relay::save2db missing staged data. Run relay.collect or relay.stage_from_db first.');
    end
    stage = getappdata(0, 'GW_RELAY_STAGE');
  end

  % Ensure output directory exists
  outdir = fileparts(filePath);
  if ~isempty(outdir) && ~exist(outdir, 'dir')
    mkdir(outdir);
  end

  % Save to file with -v7.3 nocompression
  relay_stage = stage;
  try
    save(filePath, 'relay_stage', '-v7.3', '-nocompression');
    info = struct();
    info.ok = true;
    info.message = sprintf('Successfully saved to %s', filePath);
    info.filePath = filePath;
    info.savedAt = datestr(now, 'yyyy-mm-dd HH:MM:SS');
  catch ME
    info = struct();
    info.ok = false;
    info.message = ['Failed to save: ' ME.message];
    info.filePath = filePath;
    info.error = ME;
  end

end
