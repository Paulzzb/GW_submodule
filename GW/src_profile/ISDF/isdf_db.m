%
% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
function isdf_db()
% DEMO: create a minimal ISDF database and read/write data fast

root = fullfile(pwd, "ISDF_DB");
isdf_db_init(root);


% ---- create meta skeleton (desc/xga/pmu/hVh + reserved G0set) ----
desc = struct();
desc.name = "ISDF database (skeleton)";
desc.created_at = char(datetime("now"));
desc.notes = "Fill in later";

meta = db_create(root, desc);
end


