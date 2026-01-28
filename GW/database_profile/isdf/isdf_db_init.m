% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/28 ZZ
%
function isdf_db_init(root)
% ============================================================
  if ~exist(root, "dir"); mkdir(root); end
  entriesDir = fullfile(root, "entries");
  if ~exist(entriesDir, "dir"); mkdir(entriesDir); end
  
  idxPath = fullfile(root, "index.json");
  if ~exist(idxPath, "file")
    index = struct();
    index.version = 1;
    index.created_at = char(datetime("now"));
    index.entries = struct(); % entries.(key) = struct(...)
    isdf_save_json(idxPath, index);
  end
end
