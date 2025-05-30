clear all;
close all;
%
INFO_DIR = 'Si8/';

cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/'];
KS_DIR = ['D:/kssolvGW/GW_submodule/hefeikssolv/'];
GW_DIR = ['D:/kssolvGW/GW_submodule/GW/'];
GWsetfile = [CPATH, INFO_DIR, 'setGW_fullfreq_cd.m'];
INTER_DIR = [CPATH, INFO_DIR, 'IntermediateFiles/'];
GWinput = [INTER_DIR, 'scinfo.mat'];
 

% 1. Managing paths.
cd(KS_DIR);
KSSOLV_startup;
cd(GW_DIR);
gw_startup;
cd(CPATH);

% 2. Load
load(GWinput);

% 3. run
run(GWsetfile);
options_GW.inputfile = GWinput;
[GWinput, optionsGW] = kssolv2GW(options_GW, mol);


%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    % 3.3 Read qe and ks
%    qepath = 'Si8/Si8.save'
%    [sys, extra_info] = read_qe_gw_bgw(qepath);
%    
%    
%    % use ggrid to find corresponding
%    qeinfo.ggrid = double(extra_info.mill{1});
%    ksinfo.ggrid = GWinput.coulG(:, 1:3);
%    ind_bgw2ks = zeros(size(GWinput.coulG, 1), 1);
%    for i = 1:size(GWinput.coulG, 1)
%        ind_bgw2ks(i) = findind(ksinfo.ggrid(i, :), qeinfo.ggrid);
%    end
%    
%    % change data
%    nb = optionsGW.Constant.nv + optionsGW.Constant.nc;
%    GWinput.ev = extra_info.ev(1:nb);
%    % GWinput.coulG = extra_info.coulG;
%    GWinput.Z = extra_info.X0.wavefuncell{1}.psi(ind_bgw2ks, 1:nb);
%    optionsGW.Groundstate.ind_bgw2ks = ind_bgw2ks;
%   
%   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   vcoul = readmatrix('./4-sigma2/vcoul');
%   ksinfo.ggrid = GWinput.coulG(:, 1:3);
%   bgw.ggrid = vcoul(:, 4:6);
%   ind_bgw2ks_ = zeros(size(GWinput.coulG, 1), 1);
%   for i = 1:size(GWinput.coulG, 1)
%       ind_bgw2ks_(i) = findind(ksinfo.ggrid(i, :), bgw.ggrid);
%   end
%   optionsGW.Groundstate.ind_bgw2ks_ = ind_bgw2ks_;
%   
%   wfnkq = readmatrix('./4-sigma2/wfn');
%   wfnkq = wfnkq(:, 1) + 1i * wfnkq(:, 2);
%   ng = size(GWinput.Z, 1);
%   wfnkq = reshape(wfnkq, ng, []);
%   GWinput.Z = wfnkq(ind_bgw2ks, 1:nb);
%   
%   ev = load('./4-sigma2/ev');
%   GWinput.ev = ev(1:nb) / optionsGW.Constant.ry2ev;



% 4. Save
save([INTER_DIR, 'GWinput.mat'], 'GWinput', 'optionsGW');



function ind = findind(rowvec, mat)
  [m, n] = size(mat);
  if (length(rowvec) ~= n)
    error('Size does not match!');
  end
  for i = 1:m
    if (norm(rowvec - mat(i, :)) < 1e-12)
      ind = i;
      return;
    end
  end
  ind = -1;
end



