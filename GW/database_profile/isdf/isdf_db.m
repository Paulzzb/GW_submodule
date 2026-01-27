%
% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
function isdf_id = isdf_db()
% DEMO: create a minimal ISDF database and read/write data fast

root = fullfile(pwd, "ISDF_DB");
if ~exist(root, "dir"); mkdir(root); end

% ---- create meta skeleton (desc/xga/pmu/hVh + reserved G0set) ----
desc = struct();
desc.name = "ISDF database (skeleton)";
desc.created_at = char(datetime("now"));
desc.notes = "Fill in later";

meta = db_create(root, desc);
end


