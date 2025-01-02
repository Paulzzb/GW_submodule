function Omega2 = wpeff(ksinfo, options, igcol, qindex)
%% Calculate effective plasma frequencies(squared).
%% Specifically, for a give igcol, calculate all igrow for it.
%% Omega(G,G`)^2 = wp^2 * [rho(G-G`)/rho(0)] * (q+G).(q+G`)*vc(q+G)/(8pi)
%% Units are eV^2.

%% Vars that need to change when inducing k-points
%% Extra input:
%%   isrtrq         g --> g+q, size ng * 1, use it to replace all igadd and 
%%                  igpadd.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
testflag = true;
testflag = false;

nameConstants = fieldnames(options.Constant);
for i = 1:numel(nameConstants)
    fieldname = nameConstants{i};
    value = options.Constant.(fieldname);    
    % Now, you can use the variable "fieldname" and "value" in your function
    strEval = sprintf('%s = %.16f;', fieldname, value);
    eval(strEval);
%    fprintf('Field: %s, Value: %s\n', fieldname, num2str(value));
end

if ~exist('qindex') % for situations without k-points
  qindex = 1;
end



%% Initialize, get info from ksinfo
mol = ksinfo.mol;    
coulG = ksinfo.coulG(:, 4);
coulG(1) = ksinfo.coulG0;
coulG = coulG / ksinfo.vol;
rho = ksinfo.rho;
G_index = ksinfo.coulG(:, 1:3);
idxnz = ksinfo.idxnz;
gvec = ksinfo.gvecrho;
qk = ksinfo.qk;
% rvec_index = ksinfo.rvec_index;
% gvec_index = ksinfo.gvec_index;

%% Allocate space
Omega2  = zeros(ng, 1); % Same as wpmtx in wpeff output
precalc = zeros(ng, 3);
fact_wpeff = 16*pi*ry2ev^2/vol;


%% CONSTANT
coulfact_sigma_main = 8 * pi / ksinfo.vol;      
% sigma_main line 1038, coulfact = 8D0*PI_D/(dble(gr%nf-sig%nq0+1)*crys%celvol)
qnorm = dot(qk(qindex, :), qk(qindex, :)); % for k-point situation, modify here.
plasma_omega = sqrt(4 * pi * rho(1)) * 2; % 先认为某种原因下，这里确实应该是16
% grid = Ggrid(mol);
% idxnz = grid.idxnz;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start calculation

%% Check q_is_not_zero
if qnorm > TOL_ZERO
  q_is_not_zero = true;
else
  q_is_not_zero = false;
end

%% Loop over g and tabulate 
%% line 114 to line 133 of wpeff.f90


% precalc_bgw = readfile('precalc', '%f');
% precalc_bgw = reshape(precalc_bgw{1}, 3, []);
% precalc_bgw = precalc_bgw';
% precalc_bgw = precalc_bgw(ksinfo.bgw2ks, :);

if ~q_is_not_zero % Currently, always goes into this loop
  for igrow = 2:ng % line 122 to 131 of wpeff.f90
    igadd = igrow; 
    qg = G_index(igadd, :) + 0.0; % add qpoints here later
    % precalc(igrow, :) = fact_wpeff * coulG(igrow) * qg / coulfact_sigma_main; 
    precalc(igrow, :) = fact_wpeff * coulG(igrow) * qg / coulfact_sigma_main; 
  end % for igrow
  for igrow = 1:1 % if g --> 0
    precalc(igrow, :) = 0.0;
  end
else % for multi-kpoints situation
  for igrow = 1:ng
    igadd = igrow; % line 118, find index of q+g.
    if igadd ~= 1
      qg = G_index(igadd, :) + 0.0; % add qpoints here later.
      % precalc(igrow, :) = fact_wpeff * coulG(igrow) * qg / coulfact_sigma_main;
      precalc(igrow, :) = fact_wpeff * coulG(igrow) * qg / coulfact_sigma_main;
    else
      precalc(igrow, :) = 0.0;
    end % if igadd ~= 1
  end
end 


if (testflag  && igcol == 1)
  tocompare = readmatrix('precalc.csv');
  diff = norm(precalc - tocompare) / norm(precalc);
  fprintf("precalc difference = %.3e\n", diff);
end

igpadd = igcol; % line 142
qgp    = G_index(igpadd, :) + qk(qindex, :);

for igrow = 1:ng
  igadd = igrow; % line 150, for corresponding relations between G to q+G.
  Omega2(igadd) = 0;
  gg = G_index(igadd, :) - G_index(igpadd, :);
%  gg1D = mill2nl(gg, mol.n1, mol.n2, mol.n3);
  
%   disp([num2str(gg), ' ',num2str(gg1D),  ' ', num2str(gg1Dstar)])
  kadd = findvector(gg, gvec);
%   kadd = find(idxnz == gg1D);
%   if isempty(kadd)
%     continue;
%   end
  if testflag
    dlmwrite('findvector_inwpeff.csv', kadd, '-append');
  end
  if kadd == 0
    if testflag
      dlmwrite('rho_g_minus_gp.csv', 0, '-append');
      dlmwrite('Omega2.csv', 0, '-append');
    end
    continue;
  end
	rho_g_minus_gp = rho(kadd); %% sum over kpoints if necessary.
  if testflag
    dlmwrite('rho_g_minus_gp.csv', rho_g_minus_gp, '-append', 'precision', 16);
  end
  if (igadd ~= 1 || q_is_not_zero)
  %  Omega2(igrow) =  (qgp * bdot * precalc(igrow, :)') * rho_g_minus_gp;     
    Omega2(igadd) =  (qgp * ksinfo.bdot * precalc(igrow, :)') * rho_g_minus_gp;
  elseif igpadd == 1
    Omega2(igadd) = fact_wpeff * rho_g_minus_gp;
  else
    Omega2(igadd) = 0;
  end
  if testflag
    dlmwrite('Omega2.csv', Omega2(igadd), '-append', 'precision', 16);
  end
end % for igrow

% if ~q_is_not_zero % Currently, always false
%   for igrow = 1:1
%     igadd = igrow; % line 150, for corresponding relations between G to q+G.
%     gg = G_index(igadd, :) - G_index(igpadd, :);
%     gg1D = mill2nl(gg, mol.n1, mol.n2, mol.n3);
%     kadd = find(idxnz == gg1D);
%     if ~isempty(kadd)
%       rho_g_minus_gp = rho(kadd); %% sum over kpoints if necessary.
%     %  Omega2(igrow) =  (qgp * bdot * precalc(igrow, :)') * rho_g_minus_gp;     
%       Omega2(igrow) =  fact_wpeff * rho_g_minus_gp;
%     end
% %    rho_g_minus_gp = rho(rvec_index(gg1D)); %% sum over kpoints if necessary.
% %    Omega2(igrow) = 0; % line 185, wpeff
%   end     
%   for igrow = 2:ng
%     igadd = igrow; % line 150, for corresponding relations between G to q+G.
%     gg = G_index(igadd, :) - G_index(igpadd, :);
%     gg1D = mill2nl(gg, mol.n1, mol.n2, mol.n3);
%     
% %    disp([num2str(gg), ' ',num2str(gg1D),  ' ', num2str(gg1Dstar)])
% %    kadd = findvector(gg, gvec, options.setting_method);
%     kadd = find(idxnz == gg1D);
%     if ~isempty(kadd)
%       rho_g_minus_gp = rho(kadd); %% sum over kpoints if necessary.
%     %  Omega2(igrow) =  (qgp * bdot * precalc(igrow, :)') * rho_g_minus_gp;     
%       Omega2(igrow) =  (qgp * ksinfo.bdot * precalc(igrow, :)') * rho_g_minus_gp;
%     end     
%   end % for igrow
% else % if q_is_not_zero
%   for igrow = 2:ng
%     igadd = igrow; % line 150, for corresponding relations between G to q+G.
%     gg = G_index(igadd, :) - G_index(igpadd, :);
%     gg1D = mill2nl(gg, mol.n1, mol.n2, mol.n3);
%     
% %    disp([num2str(gg), ' ',num2str(gg1D),  ' ', num2str(gg1Dstar)])
% %    kadd = findvector(gg, gvec, options.setting_method);
%     kadd = find(idxnz == gg1D);
%     if ~isempty(kadd)
%       rho_g_minus_gp = rho(kadd); %% sum over kpoints if necessary.
%     %  Omega2(igrow) =  (qgp * bdot * precalc(igrow, :)') * rho_g_minus_gp;     
%       Omega2(igrow) =  (qgp * ksinfo.bdot * precalc(igrow, :)') * rho_g_minus_gp;
%     end     
%   end % for igrow
%  for igrow = 1:ng
%    igadd = igrow; % line 150, for corresponding relations between G to q+G.
%    gg = G_index(igadd, :) - G_index(igpadd, :);
%    gg1D = ((gg(1)+get(mol,'n1')/2)*get(mol,'n2') ...
%           + gg(2)+get(mol,'n2')/2)*get(mol,'n3') + gg(3) + get(mol,'n3')/2+1;
%    rho_g_minus_gp = rho(rvec_index(gg1D)); %% sum over kpoints if necessary.
%    if igadd ~= 1 % line 166
%      Omega2(igrow) = (qgp * precalc(igrow, :)') * rho_g_minus_gp;     
%    else % line 177. %%%%%% HERE MAY STILL WRONG!!!!
%      Omega2(igrow) = fact_wpeff * rho_g_minus_gp;
%    end
%  end % for igrow
% end % if q_is_not_zero

return % function Omega = wpeff()
end % function Omega = wpeff()

% function nl=mill2nl(mill,n1,n2,n3)
% %Convert mill index to nl (the index of the full G array)
% assert(size(mill,2)==3,'Sencond dimension of mill should be 3!')
% m1 = mill(:,1); m2 = mill(:,2); m3 = mill(:,3);
% m1= m1 + ((m1<0) * n1);
% m2= m2 + ((m2<0) * n2);
% m3= m3 + ((m3<0) * n3);
% nl= m1 + m2 * n1 + m3 * n1 * n2 + 1;
% end

% function iout = findvector(kk, gvec)
% 	iout = mill2nl(kk, gvec.fftgrid(1), gvec.fftgrid(2), gvec.fftgrid(3));
% 	if (iout >= 1 && iout <= gvec.nfftgridpts)
%     iout = gvec.index_vec(iout);
% 		if iout >= 1 && iout <= gvec.ng
% 			if (any(kk ~= gvec.components(iout, :)))
% 				iout = 0;
% 			end
% 		else
% 				iout = 0;
% 		end
% 	else
% 		iout = 0;
% 	end
%   
% 	return
% end
