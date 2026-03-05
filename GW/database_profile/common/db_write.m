% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
%
% Last modified: 2026/01/28 ZZ
function meta = db_write(root, meta, ID, key, A)
% ================================================================
% Fast write array A to data/<key>.bin and record dtype/shape in meta.json
% Supports real/complex for single/double (extend if you need more).

sID = num2str(ID);
dbfull = ["data"+ sID];
dataDir = fullfile(root, dbfull);
if ~exist(dataDir, "dir"); mkdir(dataDir); end

info = struct();
info.key = key;
info.shape = size(A);
info.isComplex = ~isreal(A);

cls = class(A);
info.class = cls;

binPath = fullfile(dataDir, key + ".bin");

fid = fopen(binPath, "Wb");
assert(fid > 0, "Cannot open %s for writing.", binPath);
cleanObj = onCleanup(@() fclose(fid));

if info.isComplex
  % store as interleaved [real; imag] blocks (contiguous)
  % layout: real(A(:)) then imag(A(:))
  fwrite(fid, real(A(:)), cls);
  fwrite(fid, imag(A(:)), cls);
else
  fwrite(fid, A(:), cls);
end

info.file = "data" + sID + "/" + key + ".bin";
meta.datasets.(dbfull).fields.(key) = info;

db_save(root, meta);

end

