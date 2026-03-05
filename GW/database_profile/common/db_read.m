% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
%
% Last modified: 2026/01/28 ZZ
function A = db_read(root, ID, meta, key)
% ================================================================
% Fast read array from the recorded file and reshape
sID = num2str(ID);
dataID = "data"+sID;
assert(isfield(meta.datasets.(dataID).fields, key), "Key '%s' not found in meta.", key);
info = meta.datasets.(dataID).fields.(key);
assert(isfield(info, "file"), "Key '%s' exists but has no file stored.", key);

binPath = fullfile(root,char(info.file));

fid = fopen(binPath, "rb");
assert(fid > 0, "Cannot open %s for reading.", binPath);
cleanObj = onCleanup(@() fclose(fid));

cls = char(info.class);
n = prod(info.shape);

if isfield(info, "isComplex") && info.isComplex
    re = fread(fid, n, cls);
    im = fread(fid, n, cls);
    A = reshape(complex(re, im), info.shape);
else
    x = fread(fid, n, cls);
    A = reshape(x, info.shape);
end
end %end func
