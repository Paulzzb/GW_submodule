%
% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
function A = db_read(root, meta, key)
% ================================================================
% Fast read array from the recorded file and reshape

assert(isfield(meta.fields, key), "Key '%s' not found in meta.", key);
info = meta.fields.(key);
assert(isfield(info, "file"), "Key '%s' exists but has no file stored.", key);

binPath = fullfile(root, char(info.file));

fid = fopen(binPath, "Rb");
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
en %end func
