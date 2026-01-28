% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
%
% Last modified: 2026/01/28 ZZ
function meta = db_create(root, desc, Ndb)
% Create DB with Ndb datasets: data1, data2, data3, ...data($Ndb)
% Each dataset has the same field placeholders: xga/pmu/hVh/G0set
if ~exist(root, "dir"); mkdir(root); end

meta = struct();
meta.version = 2;
meta.desc = desc;

% datasets container
meta.datasets = struct();

% create three dataset dirs and placeholders
for k = 1:Ndb
  dname = sprintf("data%d", k);
  dataDir = fullfile(root, dname);
  if ~exist(dataDir, "dir"); mkdir(dataDir); end

  ds = struct();
  ds.path = dname; % relative path
  ds.fields = struct();
  ds.fields.xga = struct();    % placeholder
  ds.fields.pmu = struct();    % placeholder
  ds.fields.hVh = struct();    % placeholder
  ds.fields.G0set = struct();  % placeholder (reserved)

  meta.datasets.(dname) = ds;
end

db_save(root, meta);
end
