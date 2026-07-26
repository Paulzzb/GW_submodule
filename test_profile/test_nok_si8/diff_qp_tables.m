function outPath = diff_qp_tables(fileA, fileB, outFile)
%DIFF_QP_TABLES  Difference two qp_*.dat tables (full-frequency two-line layout).
%
%   outPath = diff_qp_tables(fileA, fileB)
%   outPath = diff_qp_tables(fileA, fileB, outFile)
%
% Parses the format written by qp_cohsex_fout (frequency_dependence == 2):
%   line1: n, Emf, Eo, X, Re SX-X, Re CH, Re Sig, Vxc, Re Eqp0
%   line2: (pad) Im SX-X, Im CH, Im Sig, (pad) Im Eqp0
%
% Output columns are (fileB - fileA). Defaults in this folder:
%   fileA = qp_cohsex.dat, fileB = qp_ff_dir.dat,
%   outFile = qp_ff_dir_minus_cohsex.dat

if nargin < 1 || isempty(fileA)
  fileA = fullfile(fileparts(mfilename('fullpath')), 'qp_cohsex.dat');
end
if nargin < 2 || isempty(fileB)
  fileB = fullfile(fileparts(mfilename('fullpath')), 'qp_ff_dir.dat');
end
if nargin < 3 || isempty(outFile)
  outFile = fullfile(fileparts(mfilename('fullpath')), 'qp_ff_dir_minus_cohsex.dat');
end

[ReA, ImA] = local_read_qp_ff_table(fileA);
[ReB, ImB] = local_read_qp_ff_table(fileB);
if size(ReA, 1) ~= size(ReB, 1)
  error('diff_qp_tables:nband', 'Row count mismatch: %s has %d bands, %s has %d.', ...
    fileA, size(ReA, 1), fileB, size(ReB, 1));
end
dRe = ReB - ReA;
dIm = ImB - ImA;
dRe(:, 1) = ReB(:, 1);

fid = fopen(outFile, 'w');
if fid == -1
  error('diff_qp_tables:open', 'Cannot write: %s', outFile);
end
cln = onCleanup(@() fclose(fid));
fprintf(fid, '# delta = fileB - fileA\n');
fprintf(fid, '# fileA = %s\n', fileA);
fprintf(fid, '# fileB = %s\n', fileB);
fprintf(fid, '   n      dEmf       dEo         dX    dRe SX-X    dRe CH   dRe Sig      dVxc   dRe Eqp0\n');
fprintf(fid, '                               dIm SX-X    dIm CH  dIm Sig            dIm Eqp0\n');
for i = 1:size(dRe, 1)
  fprintf(fid, '%4d%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f%12.6f \n', dRe(i, :));
  fprintf(fid, '%40s%12.6f%12.6f%12.6f            %12.6f \n', '', dIm(i, :));
end
if nargout > 0
  outPath = outFile;
end
fprintf('Wrote %s\n', outFile);
end

function [ReBlock, ImBlock] = local_read_qp_ff_table(path)
lines = local_read_nonempty_lines(path);
if numel(lines) < 3
  error('diff_qp_tables:empty', 'Too few lines in %s', path);
end
i0 = 1;
while i0 <= numel(lines) && startsWith(strtrim(lines{i0}), '#')
  i0 = i0 + 1;
end
hdr = 2;
i0 = i0 + hdr;
body = lines(i0:end);
npair = floor(numel(body) / 2);
ReBlock = zeros(npair, 9);
ImBlock = zeros(npair, 4);
for p = 1:npair
  L1 = body{2 * p - 1};
  L2 = body{2 * p};
  r = sscanf(L1, '%d %f %f %f %f %f %f %f %f');
  if numel(r) ~= 9
    error('diff_qp_tables:parseRe', 'Line %d: expected 9 numbers (n + 8 Re), got %d in:\n%s', ...
      2 * p - 1, numel(r), L1);
  end
  ReBlock(p, :) = r(:).';
  i2 = sscanf(L2, '%f %f %f %f');
  if numel(i2) ~= 4
    error('diff_qp_tables:parseIm', 'Line %d: expected 4 Im numbers, got %d in:\n%s', ...
      2 * p, numel(i2), L2);
  end
  ImBlock(p, :) = i2(:).';
end
end

function lines = local_read_nonempty_lines(path)
fid = fopen(path, 'r');
if fid == -1
  error('diff_qp_tables:read', 'Cannot open: %s', path);
end
cln = onCleanup(@() fclose(fid));
raw = textscan(fid, '%s', 'Delimiter', '\n', 'Whitespace', '');
raw = raw{1};
lines = cell(0, 1);
for k = 1:numel(raw)
  s = strtrim(raw{k});
  if isempty(s)
    continue;
  end
  lines{end + 1, 1} = raw{k}; %#ok<AGROW>
end
end
