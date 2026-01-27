% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): (AUTHOR LIST) 
%
function meta = db_create(root, desc)
% ================================================================
% Create database folder structure + meta.json

if ~exist(root, "dir"); mkdir(root); end
dataDir = fullfile(root, "data");
if ~exist(dataDir, "dir"); mkdir(dataDir); end

meta = struct();
meta.version = 1;
meta.desc = desc;

% required fields
meta.fields = struct();
meta.fields.xga = struct();   % placeholder
meta.fields.pmu = struct();   % placeholder
meta.fields.hVh = struct();   % placeholder

% reserved fields
meta.fields.G0set = struct(); % placeholder (reserved)

db_save(root, meta);
end
