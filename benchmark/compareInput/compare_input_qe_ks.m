function compare_input_qe_ks(kspath, qepath)
% This function compare the input of qe and kssolv

wavefuncpass = 1;
vxcpass = 1;
evpass = 1;
vcoulpass = 1;


% Read qe and ks
[sys, extra_info] = read_qe_gw_bgw(qepath);
load([kspath,'/GWinput.mat']);

% Now, in GWinput are kssolv gs result
% and sys/extra_info are qe gs result

qeinfo = [];
ksinfo = [];

% For k-points situation, later we fix it.

qeinfo.ggrid = double(extra_info.mill{1});
ksinfo.ggrid = GWinput.coulG(:, 1:3);
qeinfo.wf = extra_info.X0.wavefuncell{1}.psi;
ksinfo.wf = GWinput.Z;
qeinfo.ev = extra_info.ev;
ksinfo.ev = GWinput.ev;
xml=[qepath,'/vxc.dat'];
fid = dlmread(xml);
[vxcrow,vxccol]=size(fid);
step=fid(1,4);
vxc.kpoints=[];
vxc.value =[];
for i =1:step+1:vxcrow-step
vxc.kpoints=[vxc.kpoints;fid(i,1),fid(i,2),fid(i,3)];
a=[];
  for j=i+1: 1 :i+step
      a=[a,fid(j,3)];
  end
vxc.value=[vxc.value;a];
end
vxc.value=vxc.value.';
qeinfo.vxc = vxc.value;
ksinfo.vxc = GWinput.Vxc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Effective band number
nb = min(length(qeinfo.vxc), length(ksinfo.vxc));

%% indices conversion

% ks_wf(ind_ks2bgw) = bgw_wf ( hope so ... )
ind_ks2bgw = zeros(size(GWinput.coulG, 1), 1);
for i = 1:size(GWinput.coulG, 1)
    ind_ks2bgw(i) = findind(qeinfo.ggrid(i, :), ksinfo.ggrid);
end
ind_bgw2ks = zeros(size(GWinput.coulG, 1), 1);
for i = 1:size(GWinput.coulG, 1)
    ind_bgw2ks(i) = findind(ksinfo.ggrid(i, :), qeinfo.ggrid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare
% 
% ggrid
"grid dif"
norm(ksinfo.ggrid(ind_ks2bgw, :) - qeinfo.ggrid)
norm(qeinfo.ggrid(ind_bgw2ks, :) - ksinfo.ggrid)

% Compare wavefunctions
"wf diff"
norm(ksinfo.wf(ind_ks2bgw, 1:nb) - qeinfo.wf(:, 1:nb))
norm(abs(ksinfo.wf(ind_ks2bgw, 1:nb)).^2 - abs(qeinfo.wf(:, 1:nb)).^2)

% Compare Vxc
"vxc differnet"
norm(ksinfo.vxc(1:nb) - qeinfo.vxc(1:nb))




% Compare ev
"ev diff"
norm(ksinfo.ev(1:nb) - qeinfo.ev(1:nb))

% Compare vcoul
"vcoul diff"
compare_vcoul()

end


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