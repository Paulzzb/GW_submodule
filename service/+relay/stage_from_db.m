% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/23

function stage = stage_from_db(filePath)
% Load relay stage data from a DB file into temporary relay cache.
% Current fallback behavior accepts MAT files with variable relay_stage/stage.

  if nargin < 1 || isempty(filePath)
    error('relay::stage_from_db requires filePath.');
  end

  if ~exist(filePath, 'file')
    error('relay::stage_from_db cannot find file: %s', filePath);
  end

  payload = load(filePath);
  stage = [];

  if isfield(payload, 'relay_stage')
    stage = payload.relay_stage;
  elseif isfield(payload, 'stage')
    stage = payload.stage;
  else
    names = fieldnames(payload);
    for i = 1:numel(names)
      candidate = payload.(names{i});
      if isstruct(candidate)
        stage = candidate;
        break
      end
    end
  end

  if isempty(stage)
    error('relay::stage_from_db did not find valid staged data in %s', filePath);
  end

  setappdata(0, 'GW_RELAY_STAGE', stage);
end
