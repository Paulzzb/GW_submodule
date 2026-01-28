% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/28 ZZ
function meta = db_read_meta(root)
% Only load meta.json, do NOT touch any binary data

metaPath = fullfile(root, "meta.json");
assert(isfile(metaPath), "meta.json not found in %s", root);

txt = fileread(metaPath);
meta = jsondecode(txt);
end
