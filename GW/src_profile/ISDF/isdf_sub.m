% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
% Last modified: 2026/01/29 ZZ
function [ind_mu, hVh, pga] = isdf_sub(type, indicater, dbroot, varargin)
%====================================================================
% This function looks like old isdf_main function
% If indicater = 0, calculate result and save in dbroot
%              = 1, read result from dbroot and output
%====================================================================

%
switch type
  case 'vc'
    IID = 1; SID = "1";
  case 'vs'
    IID = 2; SID = "2";
  case 'ss'
    IID = 3; SID = "3";
  otherwise
    msg = sprintf("ISDF 'type' should be 'vc', 'vs', or 'ss'");
    QPerror(msg);    
end
%
meta = db_read_meta(dbroot);
%
% -------------------------------------------------------------------
% Read from database and output if indicater = 1
%
if indicater(1) == 1
  %
  msg = sprintf('Loading from %s...', dbroot);
  QPlog(msg);
  %
  ind_mu  = db_read(dbroot, IID, meta, "ind_mu");
  hVh     = db_read(dbroot, IID, meta, "hVh");
  pga = db_read(dbroot, IID, meta, "pga");
  return
end
ind_mu = 1;
hVh = 1;
pga = 1;
% -------------------------------------------------------------------
% Else, do ISDF calculation
msg = sprintf('Starting ISDF calculation...');
QPlog(msg);
%
default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end
%
if ~isempty(varargin)
  if numel(varargin) >= 1, GWinfo = varargin{1}; end
  if numel(varargin) >= 2, config = varargin{2}; end
end
%
if numel(varargin) < 2
  msg = sprintf("When implementing ISDF calculation, varargin = 6 is necessary.");
  QPerror(msg)
end
%
is_helper = false;
if config.ISDF.is_helper > 0
  warning('isdf_sub:ZetaOptional', ...
    ['zeta_mu is generally not recommended to output. ', ...
     'If you really need it, ensure cache is big enough']);
  is_helper = true;
end
% Extract data from GWinfo/config
psir = GWinfo.psir;
optionsISDF = config.ISDFCauchy;
nv = find(GWinfo.occupation > 1 - TOL_SMALL, 1, 'last');
nsum = config.SYSTEM.number_bands_in_summation;
nbmin = config.SYSTEM.energy_band_index_min;
nbmax = config.SYSTEM.energy_band_index_max;
gvec = GWinfo.gvec;
vol = GWinfo.vol;
ng = GWinfo.gvec.ng;
Dcoul = spdiags(GWinfo.coulG, 0, ng, ng);
Dcoul(1,1) = GWinfo.coulG0;
Dcoul = Dcoul * ry2ev;
%
switch IID
  case 1
    nlist = 1:nv; mlist = nv+1:nsum;
    kisdf = config.ISDF.isdf_ratio_type1;
  case 2
    nlist = 1:nv; mlist = nbmin:nbmax;
    kisdf = config.ISDF.isdf_ratio_type2;
  case 3
    nlist = 1:nsum; mlist = nbmin:nbmax;
    kisdf = config.ISDF.isdf_ratio_type3;
end
Nisdf = ceil(kisdf*sqrt(length(nlist)*length(mlist)));
optionsISDF.isdfoptions.rank = Nisdf;
tmp = "desc_type"+SID;
meta.desc.(tmp).add("Nisdf", Nisdf);
meta.desc.(tmp).add("nlist", [nlist(1), nlist(end)]);
meta.desc.(tmp).add("mlist", [mlist(1), mlist(end)]);
db_save(dbroot, meta);
%
% ===================================================================
psi = conj(psir(:, nlist));
phi = psir(:, mlist);
% Calculate this in isdf_main_
% Step 1: Compute interpolation points
msg = sprintf('Generating interpolation points...');
QPlog(msg, 2);
ind_mu = isdf_indices(psi, phi, optionsISDF);
meta = db_write(dbroot, meta, IID, "ind_xga", ind_mu);

% Step 2: Compute helper function if necessary, then hVh
hVh = zeros(Nisdf, Nisdf, 1);
if is_helper
  msg = sprintf('Constructing helper functions...');
  QPlog(msg);
  pga = isdf_kernelg(psi, phi, ind_mu, gvec, vol);
  msg = sprintf('Helper functions constructed successfully.');
  QPlog(msg);
  meta = db_write(dbroot, meta, IID, "pga", pga);
  % Step 2.2: Compute hVh
  msg = sprintf('Constructing helper V helper...');
  QPlog(msg);
  hVh = isdf_helper2hVh(pga, Dcoul, vol);
  msg = sprintf('Helper V helper constructed successfully.');
  QPlog(msg);
else
  % Step 2: Compute hVh
  msg = sprintf('Constructing helper V helper...');
  QPlog(msg);
  hVh = isdf_ind2hVh(psi, phi, ind_mu, Dcoul, gvec, vol);
  msg = sprintf('Helper V helper constructed successfully.');
  QPlog(msg);
end

 
% ===================================================================
% Save into database
msg = sprintf('Save into database...');
QPlog(msg);
meta = db_write(dbroot, meta, IID, "hVh", hVh);
msg = sprintf('Save successfully.');
QPlog(msg);

% ===================================================================
% Verification
msg = sprintf('ISDF type %1d validation\n', IID);
QPlog(msg);

switch IID
  case 1
    what = "coll";
  case 2
    if config.ISDF.is_helper
      what = "rho hf coll";
    else
      what = "hf coll";
    end
  case 3
    if config.ISDF.is_helper
      what = "rho hf coll";
    else
      what = "hf coll";
    end
end
isdf_validation(what, type, dbroot, GWinfo, config);

% ===================================================================
% Final output 
msg = sprintf('ISDF sub routine finished.');
QPlog(msg, 0);

end % function
