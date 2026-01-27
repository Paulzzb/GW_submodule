%
% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
function db_save(root, meta)
% Save meta.json (pretty + stable)
% ================================================================
metaPath = fullfile(root, "meta.json");
txt = jsonencode(meta, "PrettyPrint", true);
fid = fopen(metaPath, "w");
assert(fid > 0, "Cannot open meta.json for writing.");
cleanObj = onCleanup(@() fclose(fid));
fwrite(fid, txt, "char");
% ================================================================
end