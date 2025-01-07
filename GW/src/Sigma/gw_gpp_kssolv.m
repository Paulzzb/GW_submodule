function  eqp = gw_gpp(GWinfo, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is used to calculate GW quasiparticle energies under COHSEX approximation.
%
% Input:
% GWinfo: A structure that contains ground state calculation results
% and parameters of the molecules.
% Check gwsetup.m or README for details.
%
% options: 
%   isisdf:    true or false.
%   iscauchy:  use Cauchy integral or not, only applicable if isisdf = True.
%   optionscauchy: options when doing cauchy integral, check introduction of 
%                   COmegaCstar.m for details.
%   vcrank_mu: ISDF coefficient for valence orbitals and conduction orbitals.  
%   vsrank_mu: ISDF coefficient for valence orbitals and all orbitals.  
%   ssrank_mu: ISDF coefficient for all orbitals and all orbitals.  
%   (*rank_mu is only applicable if isisdf = true.)
%   optionsISDF: options when doing ISDF, check ./src/ISDF for details.
%     if not given, 'qrcp' is default.
%   fileName:  Output file names.
%   nv:        Number of occupied states calculated in gwsetup.m.
%   nc:        Number of unoccupied states calculated in gwsetup.m.
%   nv_ener:   Number of occupied states to calculate self energies.
%   nc_ener:   Number of unoccupied states to calculate self energies.
%   nv_oper:   Number of occupied states to manipulate operators.
%   nc_oper:   Number of unoccupied states to manipulate operators.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For test propose
% gvec_bgw = readmatrix('gvec_components.csv');
% ks2bgw = index_ks2bgw(GWinfo.coulG(:, 1:3), gvec_bgw(1:81, :));
% [~, bgw2ks] = sort(ks2bgw);
% if options.testflag == true;
%   GWinfo_bgw = setGWinfo(GWinfo, options, './');
% end


% Check field in options
if ~isfield(options, 'isisdf')
  error('Need provide whether use ISDF or not in options.isisdf!');
end
if ~isfield(options, 'fileName')
  options.fileName = 'GW_output.mat';
end
if options.isisdf
  if ~isfield(options, 'vcrank_ratio')
    error('Need provide ISDF coefficient vcrank_mu in options.vcrank_mu!');
  end
  if ~isfield(options, 'vsrank_ratio')
    error('Need provide ISDF coefficient vsrank_mu in options.vsrank_mu!');
  end
  if ~isfield(options, 'ssrank_ratio')
    error('Need provide ISDF coefficient ssrank_mu in options.ssrank_mu!');
  end
  if ~isfield(options, 'iscauchy')
    options.iscauchy = 1;
  end
  if options.iscauchy
    if isfield(options, 'optionscauchy')
      COmegaCOptions = options.optionscauchy;
    else
      COmegaCOptions.froErr = 1e-5; COmegaCOptions.MaxIter = 10;
    end
  end  
  if isfield(options, 'optionsISDF')
    optionsISDF = option.optionsISDF;
  else
    optionsISDF.isdfoptions.seed = 0;
    optionsISDF.exxmethod = 'qrcp';
  end 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import global variables
global nv nv_ener nv_oper nc nc_ener nc_oper n_oper n_ener ...
       ng nr ne vol Fouriercell ...
       TOL_ZERO TOL_SMALL INF...
       ry2ev 

% Initialization
Z     = GWinfo.Z;
ev    = GWinfo.ev;
Vxc   = GWinfo.Vxc;
vcoul = GWinfo.coulG(:,4);
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
Dcoul(1,1) = GWinfo.coulG0;
fprintf('nr = %d, ng = %d, n_oper = %d\n', nr, ng, n_oper);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Some constant, some initializing.
%% 
ZERO = 0.0;
CZERO = 0.0 + sqrt(-1) * 0.0;
%% This part of the code is useless in matlab.
%% However, we write here to remind the developers the right types of vars. 
Omega2   = CZERO * zeros(ng, 1);
wtilde_matrix = CZERO * zeros(ng, ng);
wtilde2  = CZERO;
wx_array = CZERO * zeros(1, 1); % currently, we induce no omega part, which let wx_array be size 1.
wtilde   = CZERO;
wtilde2  = CZERO;
cden     = CZERO; 
delw     = CZERO;
wdiff    = CZERO;
sch      = CZERO;
ssx      = CZERO;
ssxt     = CZERO;
scht     = CZERO;
schtt    = CZERO;
ssx_array= CZERO * zeros(1,1);
sch_array= CZERO * zeros(1,1);
asxt     = CZERO * zeros(1,1);
acht     = CZERO * zeros(1,1);
asxtemp  = CZERO * zeros(1,1);
achtemp  = CZERO * zeros(1,1);
rden     = ZERO;
delwr    = ZERO;
delw2    = ZERO;
wdiffr   = ZERO;
wxt      = ZERO;
occ      = ZERO; % this should be 0, 0.5, or 1 in real calculation, according to mtxel_cor
                 % when deciding flagocc.
sexcut   = 4.0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% New input compared with cohsex.

TOL_SMALL = 1e-6;
limitone =  1.0 / (4.0 * TOL_SMALL);
limittwo = 0.50.^2;
if isfield(options, 'gpp_brodening')
  limittwo = gpp_brodening .^ 2;
end
isCPLX = true;% ifdef CPLX

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Normalize the wavefunc in Fourier space.
startZ = tic;
for iv = 1:nv+nc
  Z(:,iv) = Z(:,iv) / (norm(Z(:,iv)));
end

timeforZ = toc(startZ)
startC2R = tic;
psir = zeros(nr, nv+nc);
% if options.testflag == true
%   Z = Z(GWinfo.BGW2ks, :);
% end
psir = F' * (Z(:, 1:nv+nc));
if options.testflag == true
  Z = Z(GWinfo.ks2BGW, :);
end
timeforC2R = toc(startC2R);

% Start our GW calculation
if (options.isisdf)
  if isfield(options, 'vcrank_ratio') & isfield(options, 'vcrank_ratio') ...
     & isfield(options, 'vcrank_ratio')
    vcrank_mu = ceil(sqrt(nv_oper*nc_oper) * options.vcrank_ratio);
    vsrank_mu = ceil(sqrt(nv_oper*n_ener)  * options.vsrank_ratio);
    ssrank_mu = ceil(sqrt(n_ener *n_ener)   * options.ssrank_ratio);
  else
    error('ISDF coefficient not provided, check options again!');
  end
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Start ISDF
  % Generally, we have relations 
  %     Mr{xx} = xxrzeta_mu * Mrxx_mu;
  startISDF = clock;
  vcrzeta_mu = zeros(nr, vcrank_mu);
  vcgzeta_mu = zeros(ng, vcrank_mu);
  
  tic
  optionsISDF.isdfoptions.rank = vcrank_mu;
  [vcrzeta_mu, vcind_mu] = isdf(sqrt(vol) * conj(psir(:,nv-nv_oper+1:nv)), (psir(:,nv+1:nv+nc_oper)), optionsISDF);
  %[vcrzeta_mu, vcind_mu] = isdf_gw(psir(:,nv-nv_oper+1:nv), conj(psir(:,nv+1:n_oper)), vcrank_mu);
  vcgzeta_mu = F * vcrzeta_mu;
  timeforVC = toc
  clear vcrzeta_mu;
  
  vsrzeta_mu = zeros(nr, vsrank_mu);
  vsgzeta_mu = zeros(ng, vsrank_mu);
  tic;
  optionsISDF.isdfoptions.rank = vsrank_mu;
  [vsrzeta_mu, vsind_mu] = isdf(sqrt(vol) * conj(psir(:,nv-nv_oper+1:nv)), (psir(:,nv-nv_ener+1:nv+nc_ener)), optionsISDF);
  % [vsrzeta_mu, vsind_mu] = isdf_gw(psir(:,nv-nv_oper+1:nv), conj(psir(:,nv-nv_oper+1:n_oper)), vsrank_mu);
  vsgzeta_mu = F * vsrzeta_mu;
  timeforVS = toc
  clear vsrzeta_mu;
  
  ssrzeta_mu = zeros(nr, vsrank_mu);
  ssgzeta_mu = zeros(ng, ssrank_mu);
  tic;
  optionsISDF.isdfoptions.rank = ssrank_mu;
  [ssrzeta_mu, ssind_mu] = isdf(sqrt(vol) * conj(psir(:,nv-nv_oper+1:nv+nc_oper)), (psir(:,nv-nv_ener+1:nv+nv_ener)), optionsISDF);
  % [ssrzeta_mu, ssind_mu] = isdf_gw(psir(:,nv-nv_oper+1:n_oper), conj(psir(:,nv-nv_oper+1:n_oper)), ssrank_mu);
  ssgzeta_mu = F * ssrzeta_mu;
  timeforSS = toc
  timeforISDF = etime(clock, startISDF);
  clear ssrzeta_mu;
  % ISDF DONE

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Manipulate operators

  % Calculate C*Omega^(-1)*C', use direct method or cauchy integral.
  startCOmegaC = tic; 

  if options.iscauchy
    Phi = sqrt(vol) * psir(vcind_mu, nv-nv_oper+1:nv); 
    Psi = psir(vcind_mu, nv+1:nc_oper+nv); 
    COmegaCOptions.froErr = 1e-5; COmegaCOptions.MaxIter = 10;
    evOcc = ev(nv-nv_oper+1:nv);
    evUnocc = ev(nv+1:nv+nc_oper);
    
    [final_results, time_c, relError] = ...
          COmegaCstar(Phi, Psi, evOcc, evUnocc, COmegaCOptions);
    clear Phi Psi evOcc evUnocc;
  else
    Mrvc_mu = prod_states_gw_(sqrt(vol) * psir(vcind_mu,nv-nv_oper+1:nv),...
                   (psir(vcind_mu,nv+1:nv+nc_oper)));
    Eden = (kron(ev(nv-nv_oper+1:nv), ones(nc_oper,1)) ... 
                - kron(ones(nv_oper, 1), ev(nv+1:nv+nc_oper)));
    final_results = Mrvc_mu * diag(Eden) * Mrvc_mu'; 
    clear Mrvc_mu Eden;
  end
  timeforCOmegaC = toc(startCOmegaC);
  

  % Calculate ( (C*Omega*C)^(-1) ...
  %                               - vcgzeta_mu' * Dcoul * vcgzeta_mu)^(-1) 
  invepsg_main = inv( inv(final_results) ./ 4 - vcgzeta_mu' * Dcoul * vcgzeta_mu );

  % Calculate operator \Sigma_{COH}, \Sigma_{SEX_X}, and \Sigma_{X} formally.
  Phivs = sqrt(vol) * psir(vsind_mu, nv-nv_oper+1:nv); 
  Phiss = sqrt(vol) * psir(ssind_mu, nv-nv_oper+1:nv+nc_oper); 
  vcDcoulvs = vcgzeta_mu' * Dcoul * vsgzeta_mu;
  vcDcoulss = vcgzeta_mu' * Dcoul * ssgzeta_mu;
  W1_mu = vcDcoulvs' * (invepsg_main * vcDcoulvs);
  W1_mu_1 = vcDcoulss' * (invepsg_main * vcDcoulss);
  Sigma_sex_x = W1_mu .* (Phivs * Phivs'); 
  Sigma_coh   = W1_mu_1 .* (Phiss * Phiss');
  clear W1_mu W1_mu_1 vcDcoulvs vcDcoulss;
  vsDcoulvs = vsgzeta_mu' * Dcoul * vsgzeta_mu;
  Sigma_x     = vsDcoulvs .* (Phivs * Phivs'); 
  clear vsDcoulvs;


  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate self-energies
  Psivs = psir(vsind_mu, nv-nv_ener+1:nv+nc_ener); 
  Psiss = psir(ssind_mu, nv-nv_ener+1:nv+nc_ener); 
  
  Esx_x = Psivs' * Sigma_sex_x * Psivs;
  Ex    = Psivs' * Sigma_x     * Psivs;
  Ech   = Psiss' * Sigma_coh   * Psiss;
  clear Psivs Psiss;


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate energies
  GWenergy.ev = ev(nv-nv_ener+1:nv+nc_ener) * ry2ev;
  GWenergy.Ex = real(diag(0.0 - Ex) * ry2ev);
  GWenergy.Esx_x = real(diag(0.0 - Esx_x) * ry2ev);
  GWenergy.Ech = real(diag(0.5 * Ech * ry2ev));
  GWenergy.Sig = real(GWenergy.Ex + GWenergy.Esx_x + GWenergy.Ech);
  GWenergy.Vxc = Vxc(nv-nv_ener+1:nv+nc_ener) * ry2ev;
  GWenergy.eqp = GWenergy.ev - GWenergy.Vxc + GWenergy.Sig;

else % No ISDF version
  
  Mgvc = zeros(ng, nv_oper * nc_oper);
	for ind_nv = nv-nv_oper+1:nv
		ind_nv_true = ind_nv + (nv - nv_oper);
		Mgvc(:, (1:nc_oper) + (ind_nv_true-1) * (nc_oper)) ...
		 = GWinfo.aqs{ind_nv}(:, nv+1:nv+nc_oper);
	end
	Mgvc = conj(Mgvc);
  fprintf('F-norm of Mgvc = %f.\n', norm(Mgvc, 'fro'));
  Eden = (1) ./ (kron(ev(nv-nv_oper+1:nv), ones(nc_oper,1)) ... 
              - kron(ones(nv_oper, 1), ev(nv+1:nv+nc_oper)));
  % Indeed chi0
	% for some reason, the result is inversed ... I believe is that sqrt(-1)
  inveps = 4 * Mgvc * (diag(Eden) * Mgvc') / (GWinfo.vol);
  fprintf('F-norm of chi = %f.\n', norm(inveps, 'fro'));
  % clear Mrvc Mgvc;
  % Indeed epsilon
	Dcoul(1) = 0;
  inveps = eye(ng) - Dcoul * inveps;
  % Real inveps
  inveps = inv(inveps);
	I_inveps = eye(ng) - inveps;
	
	% Check correctness of I_inveps
%   I_eps_array = readfile('I_eps_array', '(%f, %f)');
%   I_eps_array = complex(I_eps_array{1}, I_eps_array{2});
%   I_eps_array = reshape(I_eps_array, ng, []);
%   fprintf('F-norm of I_inveps = %f.\n', norm(I_inveps, 'fro'));
%   fprintf('F-norm of I_inveps from bgw = %f.\n', norm(I_eps_array, 'fro'));
%   fprintf('F-norm difference of eps = %f.\n', norm(I_inveps(ks2bgw, ks2bgw)-I_eps_array, 'fro'));
  clear inveps;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% GPP constructions of inveps, 
  %% check in Sigma/mtxel_cor.f90, subroutine sigma_gpp for details.
  
	Dcoul(1) = GWinfo.coulG0;
  startwtilde = tic;
  wtilde_array = zeros(ng, ng) * CZERO;
  flagscalarize = false;
  if flagscalarize == true
    wtilde2_temp = zeros(ng, 1) * CZERO;
  else
    wtilde2_temp = CZERO;
  end
	qindex = 1;
% 	if (exist('wpmtx', 'file') == 2)
%     Omega2_bgw = readfile('wpmtx', '(%f, %f)');
% 		Omega2_bgw = Omega2_bgw{1} + Omega2_bgw{2} * sqrt(-1);
% 		ncouls = sqrt(length(Omega2_bgw));
% 		Omega2_bgw = reshape(Omega2_bgw, sqrt(length(Omega2_bgw)), []);
%     Omega2_bgw = Omega2_bgw(bgw2ks, bgw2ks);
% 	end
   

  for igcol = 1:ng
    Omega2 =  wpeff(GWinfo, igcol, qindex, options.setting_method); % There should not be 1/vol here...
%     if (exist('wpmtx', 'file') == 2)
% 			Omega2_diff = Omega2 - Omega2_bgw(:, igcol);
% 	    if (norm(Omega2_diff) / norm(Omega2)) >= 1e-3
% 				fprintf('Relative norm difference of Omega2 = %f, igcol = %d.\n', (norm(Omega2_diff) / norm(Omega2)), igcol);
% 			error('Omega2 calculation is not right!');
% 		  end
% 		end
    if igcol == 11
      disp('Need that seconds for 10 iterations')
      toc(startwtilde);
    end
    if flagscalarize == false
      for igrow = 1:ng
        I_epsggp = I_inveps(igrow, igcol);
        if (abs(I_epsggp) <= TOL_SMALL); continue; end; % line 651, mtxel_cor
        if (abs(Omega2(igrow)) <= TOL_SMALL); continue; end %line 666, mtxel_cor
        if isCPLX
          wtilde2_temp = Omega2(igrow) / I_epsggp;
          lambda = abs(wtilde2_temp);
          if (lambda <= TOL_SMALL); continue; end
          phi = atan2(imag(wtilde2_temp), real(wtilde2_temp));
          if (abs(cos(phi)) <= TOL_SMALL); continue; end
          wtilde2 = lambda / cos(phi);
        else % is isCPLX
          wtilde2 = Omega2(igrow) / I_epsggp;
          if (abs(wtilde2) <= TOL_SMALL); continue; end
        end % is isCPLX
        if (real(wtilde2) < 0) % line 697 of mtxel_cor
          wtilde = CZERO + 1.0 / TOL_ZERO; % currently, we use mode 0, which is default 
        else
          wtilde = CZERO + sqrt(real(wtilde2));
        end
        wtilde_array(igrow, igcol)   = wtilde;
      end % for igrow
    else % if flagscalarize
      bigOmega2 = find(abs(Omega2(:)) >= TOL_SMALL);
      bigIinveps = find(abs(I_inveps(:, igcol)) >= TOL_SMALL);
      indtoCalc = intersect(bigIinveps, bigOmega2);
      lambda = zeros(ng, 1) * CZERO;
      wtilde2_temp = zeros(ng, 1) * CZERO;
      wtilde2 = zeros(ng, 1) * CZERO;
      wtilde = zeros(ng, 1) * CZERO;
      phi = zeros(ng, 1);
      if isCPLX
        wtilde2_temp(indtoCalc) = Omega2(indtoCalc, :) ./ I_inveps(indtoCalc, igcol);
        lambda = abs(wtilde2_temp);
        bigLambda = find(lambda >= TOL_SMALL);
        indtoCalc = intersect(indtoCalc, bigLambda);
        phi(indtoCalc) = atan2(imag(wtilde2_temp(indtoCalc)), real(wtilde2_temp(indtoCalc)));
        indtoCalc = intersect(find(abs(cos(phi)) >= TOL_SMALL), indtoCalc);
        wtilde2(indtoCalc) = lambda(indtoCalc) ./ cos(phi(indtoCalc));
      else %if isCPLX
        wtilde2(indtoCalc) = Omega2(indtoCalc, :) ./ I_inveps(indtoCalc, igcol);
        indtoCalc = intersect(indtoCalc, find(abs(wtilde2) >= TOL_SMALL));
      end
      indwtilde2Lessthanzero = intersect(indtoCalc, find(real(wtilde2)<0));
      wtilde(indtoCalc) = CZERO + sqrt(real(wtilde2(indtoCalc))); 
      wtilde(indwtilde2Lessthanzero) = CZERO + 1.0 / TOL_ZERO; 
      wtilde_array(:, igcol) = wtilde;
    end
  end % for igcol
  toc(startwtilde)
%	wtilde_array_bgw = readfile('wtilde_array', '(%f, %f)');
%  wtilde_array_bgw = complex(wtilde_array_bgw{1}, wtilde_array_bgw{2});
%  wtilde_array_bgw = reshape(wtilde_array_bgw, ng, []);
%  wtilde_array_bgw = wtilde_array_bgw(bgw2ks, bgw2ks); 
%   for igcol = 1:ng
% 	  nonzero_wtilde_array = find(wtilde_array(:, igcol)); 
% 		nonzero_wtilde_array_bgw = find(wtilde_array_bgw(:, igcol)); 
% 		if length(nonzero_wtilde_array_bgw) ~= length(nonzero_wtilde_array);
% 		  fprintf('Different number of non-zero indices!\n');
% 			fprintf('%d : bgw, %d : Ours\n', length(nonzero_wtilde_array_bgw), length(nonzero_wtilde_array));
% 			warning();
% 			continue
% 		end
% 		if any(nonzero_wtilde_array ~= nonzero_wtilde_array_bgw)
% 			disp(nonzero_wtilde_array)
% 			disp(nonzero_wtilde_array_bgw)
% 			error('Even the non-zero indices are not the same in wtilde_array');
% 		else
% 			wtilde_array_diff = 1./ wtilde_array(nonzero_wtilde_array, igcol) ...
% 			                  - 1./wtilde_array_bgw(nonzero_wtilde_array_bgw, igcol);
% %		  fprintf('igcol = %d, difference in wtilde_array = %f\n', ...
% %			        igcol, norm(wtilde_array_diff));
% 			if (norm(wtilde_array_diff) / norm(1./wtilde_array(nonzero_wtilde_array, igcol))) >= 1e-5
% 				fprintf('Norm difference of 1/wtilde_array = %f.\n', ...
% 				        (norm(wtilde_array_diff) / norm(1./wtilde_array(nonzero_wtilde_array, igcol))));
% 				warning('wtilde_array seems wrong!')
% 				pause(0.25)
% 			end;
% 		end	
% 	end
  
	% transform their units to Ry based.
  % Notice that in wpeff.f90, units are ev.^2
  % wtilde2_matrix = wtilde2_matrix / (ry2ev.^2);
  % wtilde_array = wtilde_array / ry2ev;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% GPP construction of Sigma
  Esx_x = zeros(n_ener, 1);
  Ech   = zeros(n_ener, 1);
  Ex    = zeros(n_ener, 1);
  for ibandouter = nv-nv_ener+1 : nv+nc_ener % line 1167 of mtxel_cor
    ibandouter_aqsntemp = ibandouter - nv + nv_ener;
    if ibandouter == nv-nv_ener+1
      starttimeouter = tic;
    else
      toc(starttimeouter)
    end
    asxt = ZERO;
    acht = ZERO;
    if ~isfield(GWinfo, 'aqs')
		  aqsntempr = prod_states_gw_(sqrt(vol) * psir(:, ibandouter), ...
                    psir(:, nv-nv_oper+1:nv+nc_oper));
      aqsntemp = F * aqsntempr; % from mtxel.f90, 
                                % if you notice how mtxel_cor is called, 
                                % you understand that aqsm and aqsn is indeed
                                % the same...
                                % WHEN introducing k-points, you need an extra mapping
                                % to map g to g+q, which is none of my business :) 
    else
			aqsntemp = ...
			GWinfo.aqs{ibandouter_aqsntemp};
		end
		ssxflag = false;
		schflag = false;
		schttflag = false;
    asxtemp = ZERO;
    achtemp = ZERO; 
    for ibandinner = nv-nv_oper+1 : nv+nc_oper % line 1214 of mtxel_cor
%			ssx_name = ['ssx_occ', num2str(ibandouter_aqsntemp), '_', num2str(ibandinner)];
%			sch_name = ['sch_occ', num2str(ibandouter_aqsntemp), '_', num2str(ibandinner)];
%			schtt_name = ['sch_uno', num2str(ibandouter_aqsntemp), '_', num2str(ibandinner)];
%       if (exist(ssx_name, 'file')==2)
% 				ssx_bgw = readfile(ssx_name,'(%f, %f)');
% 				ssx_bgw = ssx_bgw{1} + sqrt(-1) * ssx_bgw{2};
% 				ssx_bgw = reshape(ssx_bgw, ng, ng);
% 			%	ssxflag = true;
% %			  ssx_array = zeros(size(ssx_bgw));
% 			elseif (exist('ssx_bgw')==1)
% 				clear ssx_bgw;
% 				ssxflag = false;
% 			end
%       if (exist(sch_name, 'file')==2)
% 				sch_bgw = readfile(sch_name, '(%f, %f)');
% 				sch_bgw = sch_bgw{1} + sqrt(-1) * sch_bgw{2};
% 				sch_bgw = reshape(sch_bgw, ng, ng);
% 			%	schflag = true;
% %			  sch_array = zeros(size(sch_bgw));
% 			elseif (exist('sch_bgw')==1)
% 				clear sch_bgw;
% 				schflag = false;
% 			end
%       if (exist(schtt_name, 'file')==2)
% 				schtt_bgw = readfile(schtt_name, '(%f, %f)');
% 				schtt_bgw = schtt_bgw{1} + sqrt(-1) * schtt_bgw{2};
% 				schtt_bgw = reshape(schtt_bgw, ng, ng);
% 			%	schttflag = true;
% %			  sch_array = zeros(size(sch_bgw));
% 			elseif (exist('schtt_bgw')==1)
% 				clear schtt_bgw;
% 				schttflag = false;
% 			end
			ibandinner_aqsntemp = ibandinner + nv_oper - nv;
      STATEMENT = (ibandinner <= nv); % refers to line 1221 
      if (STATEMENT) % we need adjustment of flagocc here!!!
        flagocc = true; occ = 1;
      else
        flagocc = false; occ = 0;
      end 
      achstemp = 0.0 + sqrt(-1) * 0.0;
      wx_array = -(ev(ibandinner) - ev(ibandouter_aqsntemp)) * ry2ev; % line 1452 of mtxel_cor;
      for igcol = 1:ng % line 1473, do my_igp = 1, ngpown
        ssx_array = CZERO * zeros(1,1);
        sch_array = CZERO * zeros(1,1);
        if flagocc == true
          scht = 0.0;
          ssxt = 0.0;
          wxt = wx_array;
          for igrow = 1:ng % line 1508
            wtilde = wtilde_array(igrow, igcol);
            wtilde2 = wtilde .^ 2;
            Omega2 = wtilde2 * I_inveps(igrow, igcol);
            wdiff = wxt - wtilde;
            cden = wdiff;
            rden = cden * conj(cden);
            rden = 1.0 / rden;
            delw = wtilde * conj(cden) * rden; % just help it to inverse a conplex number.
            delwr = delw * conj(delw);
            wdiffr = wdiff * conj(wdiff);
            if ((wdiffr >= limittwo)  &&  (delwr <= limitone)) % line 1536 of mtxel_cor
              sch  = delw * I_inveps(igrow, igcol);
              cden = wxt.^2 - wtilde2;
              rden = cden * conj(cden);
              rden = 1.0 / rden;
              ssx  = Omega2 * conj(cden) * rden;
            elseif (delwr >= TOL_ZERO) % ???why doing this???
              sch  = ZERO;
              cden = 4.0 * wtilde2 * (delw + 0.5);
              rden = cden * conj(cden);
              rden = 1.0 / rden;
              ssx  = -Omega2 * conj(cden) * rden * delw;
            else
              sch  = 0.0;
              ssx  = 0.0;
            end % if (wdiffr...)
            sexcutoff = sexcut * abs(I_inveps(igrow, igcol));
            if ((abs(ssx) > sexcutoff) && wxt < 0.0)
              ssx = 0.0;
            end
            ssxt = ssxt + ssx * aqsntemp(igrow, ibandinner_aqsntemp);
            scht = scht + sch * aqsntemp(igrow, ibandinner_aqsntemp);
%            if (ssxflag == true)
%				  		if abs(ssx) >= TOL_ZERO
%				  		  if (abs(ssx - ssx_bgw(igrow, igcol)) / abs(ssx) >= 1e-4)
%								  fprintf('ibandouter = %d, ibandinner = %d.\n', ibandouter, ibandinner);	
%									fprintf('ssx is not matched when igrow = %d, igcol = %d\n', igrow, igcol);
%				  				warning(num2str(abs(ssx - ssx_bgw(igrow, igcol)) / abs(ssx)));
%									pause(0.05)
%%									Demo2
%%									Demo3
%				  			end
%							elseif abs(ssx_bgw(igrow, igcol)) >= TOL_ZERO * 2
%								fprintf('ibandouter = %d, ibandinner = %d.\n', ibandouter, ibandinner);	
%								fprintf('igrow = %d, igcol = %d\n', igrow, igcol);
%							  warning('When ssx is small, from bgw is not!')
%								pause(0.05)
%%								Demo2
%%								Demo3
%				  		end
%				  	end
%            if (schflag == true)
%				  		if abs(sch) >= TOL_ZERO
%				  		  if (abs(sch - sch_bgw(igrow, igcol)) / abs(sch) >= 1e-4)
%								  fprintf('ibandouter = %d, ibandinner = %d.\n', ibandouter, ibandinner);	
%									fprintf('sch is not matched when igrow = %d, igcol = %d\n', igrow, igcol);
%				  				warning(num2str(abs(sch - sch_bgw(igrow, igcol)) / abs(sch)));
%									pause(0.05)
%%									Demo2
%%									Demo3
%				  			end
%							elseif abs(sch_bgw(igrow, igcol)) >= TOL_ZERO * 2
%								fprintf('ibandouter = %d, ibandinner = %d.\n', ibandouter, ibandinner);	
%								fprintf('igrow = %d, igcol = %d\n', igrow, igcol);
%							  warning('When sch is small, from bgw is not!')
%								pause(0.05)
%%								Demo2
%%								Demo3
%				  		end
%				  	end
					end % for igrow
				  ssx_array = ssx_array + ssxt * conj(aqsntemp(igcol, ibandinner_aqsntemp));
          sch_array = sch_array + 0.5 * scht * conj(aqsntemp(igcol, ibandinner_aqsntemp));
        else % if flagocc
          % from line 1568 to 1608
          scht = 0.0;
          ssxt = 0.0;
          wxt = wx_array;
          for igrow = 1:ng % line 1577 of mtxel_cor
            wdiff = wxt - wtilde_array(igrow, igcol);
            cden = wdiff;
            rden = cden * conj(cden);
            rden = 1.0 / rden;
            delw = wtilde_array(igrow, igcol) * conj(cden) * rden;
            delwr = delw * conj(delw);
            wdiffr = wdiff * conj(wdiff);
            schtt = delw * I_inveps(igrow, igcol) * aqsntemp(igrow, ibandinner_aqsntemp);
            if (wdiffr > limittwo  &&  delwr < limitone)
              scht = scht + schtt;
            end
%             if (schttflag == true)
% 				  		if abs(schtt) >= TOL_ZERO
% 				  		  if (abs(schtt - schtt_bgw(igrow, igcol)) / abs(schtt) >= 1e-4)
% 								  fprintf('ibandouter = %d, ibandinner = %d.\n', ibandouter, ibandinner);	
% 									fprintf('schtt is not matched when igrow = %d, igcol = %d\n', igrow, igcol);
% %									Demo4
% %									Demo5
% 				  				warning(num2str(abs(schtt - schtt_bgw(igrow, igcol)) / abs(schtt)));
% 									pause(0.05)
% 				  			end
% 							elseif abs(schtt_bgw(igrow, igcol)) >= TOL_ZERO * 2
% 								fprintf('ibandouter = %d, ibandinner = %d.\n', ibandouter, ibandinner);	
% 								fprintf('igrow = %d, igcol = %d\n', igrow, igcol);
% %								Demo4
% %								Demo5
% 							  warning('When schtt is small, from bgw is not!')
%                 pause(0.05) 
% 				  		end
% 				  	end % if (schttflag == true)

          end % for igrow = 1:ng
          sch_array = sch_array + 0.5 * scht * conj(aqsntemp(igcol, ibandinner_aqsntemp));
        end % if flagocc

        if (flagocc == true)
          asxtemp = asxtemp - ssx_array * occ * vcoul(igcol) / GWinfo.vol;
        end
        achtemp = achtemp + sch_array * vcoul(igcol) / GWinfo.vol;
      end % for igcol = 1:ng

    end % for ibandinner
    asxt = asxt + asxtemp;% line 1280, mtxel_cor
    acht = acht + achtemp;
    Esx_x(ibandouter_aqsntemp) = asxt;
    Ech(ibandouter_aqsntemp)   = acht;

    
  end % for ibandouter

  % W1 = inveps - Dcoul;
  % Esx_x = zeros(n_ener);
  % Ex     = zeros(n_ener);
  % Ech    = zeros(n_ener);
  % for ioper = nv-nv_oper+1:nv
  %   psiiopersPsi = sqrt(vol) * conj(psir(:, ioper)) .* psir(:, nv-nv_ener+1:nv+nc_ener);
  %   tmp = F * psiiopersPsi;
  %   Esx_x = Esx_x + tmp' * W1 * tmp;
  %   Ex = Ex + tmp' * Dcoul * tmp;
  %   Ech = Ech + tmp' * W1 * tmp;
  % end
  % for ioper = nv+1:nv+nc_oper
  %   psiiopersPsi = sqrt(vol) * conj(psir(:, ioper)) .* psir(:, nv-nv_ener+1:nv+nc_ener);
  %   tmp = F * psiiopersPsi;
  %   Ech = Ech + tmp' * W1 * tmp;
  % end

  % Calculate energies
  GWenergy.ev = ev(nv-nv_ener+1:nv+nc_ener) * ry2ev;
  GWenergy.Ex = real((0.0 + Ex) * ry2ev);
  GWenergy.Esx_x = real((0.0 + Esx_x) * ry2ev);
  GWenergy.Ech = real((1.0 * Ech * ry2ev));
  GWenergy.Sig = real(GWenergy.Ex + GWenergy.Esx_x + GWenergy.Ech);
  GWenergy.Vxc = Vxc(nv-nv_ener+1:nv+nc_ener) * ry2ev;
  GWenergy.eqp = GWenergy.ev - GWenergy.Vxc + GWenergy.Sig;

end


% save(['GWenergy_',num2str(vcrank_mu),'_',num2str(vsrzeta_mu),'_',num2str(ssrank_mu),'.mat'], GWenergy);
save(options.fileName, 'GWenergy') 
% save(['GWenergy', name_of_mol, '_',num2str(vcrank_mu), '_', num2str(vsrank_mu), '_', num2str(ssrank_mu), '.mat'], 'GWenergy');
% save(['timeinfo', name_of_mol], 'time*');
end


function Omega2 = wpeff(GWinfo, igcol, qindex, setting_method)
%% Calculate effective plasma frequencies(squared).
%% Specifically, for a give igcol, calculate all igrow for it.
%% Omega(G,G`)^2 = wp^2 * [rho(G-G`)/rho(0)] * (q+G).(q+G`)*vc(q+G)/(8pi)
%% Units are eV^2.

%% Vars that need to change when inducing k-points
%% Extra input:
%%   isrtrq         g --> g+q, size ng * 1, use it to replace all igadd and 
%%                  igpadd.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Import global variables
global nv nv_ener nv_oper nc nc_ener nc_oper n_oper n_ener ...
       ng ne vol Fouriercell ...
       TOL_ZERO TOL_SMALL INF...
       ry2ev ... 
       qk; 
if ~exist('qindex') % for situations without k-points
  qindex = 1;
end


%% Initialize, get info from GWinfo
coulG = GWinfo.coulG(:, 4);
coulG(1) = GWinfo.coulG0;
coulG = coulG / GWinfo.vol;
rho = GWinfo.rho;
G_index = GWinfo.coulG(:, 1:3);
idxnz = GWinfo.idxnz;
gvec = GWinfo.gvec;

% rvec_index = GWinfo.rvec_index;
% gvec_index = GWinfo.gvec_index;

%% Allocate space
Omega2  = zeros(ng, 1); % Same as wpmtx in wpeff output
precalc = zeros(ng, 3);

%% CONSTANT
fact_wpeff = 16 * pi * ry2ev.^2 * ne * (1 / ne) / vol;
coulfact_sigma_main = 8 * pi / GWinfo.vol;      
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
% precalc_bgw = precalc_bgw(GWinfo.bgw2ks, :);

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
% fprintf('Precalc difference = %f', norm(precalc_bgw - precalc, 'fro'))
% Here, THE CODE IS FOR DEBUGGING PRECALC
% if (exist('precalc', 'file') == 2)
%   precalc_bgw = readfile('precalc', '%f');
% 	if (igcol ~= precalc_bgw{1}(1))
% % 	fprintf('igcol %d not matched with igp %d.\n', igcol, precalc_bgw{1}(1));
%     ;
% 	else
% 		precalc_bgw = precalc_bgw{1}(2:end);
% 		precalc_bgw = reshape(precalc_bgw, 3, []);
% 		precalc_bgw = precalc_bgw';
%     % Check whether precalc_bgw is same with precalc
% 		precalc_diff = precalc - precalc_bgw;
% 		if norm(precalc_diff, 'fro') / norm(precalc, 'fro') >= 1e-6
% 			fprintf('relative error of precalc is %f! \n', norm(precalc_diff, 'fro') / norm(precalc, 'fro'));
% 			error('precalc is not the same!');
% 		end
% 		clear precalc_diff precalc_bgw
% 	end
% end


% tmp=1:ng;
% G_tmp = G(repmat(tmp',1,ng),:) - G(repmat(tmp,ng,1),:);
% clear tmp;
% corr = reshape(((g_tmp(:,1)+get(mol,'n1')/2)*get(mol,'n2') ...
%        + g_tmp(:,2)+get(mol,'n2')/2)*get(mol,'n3') + g_tmp(:,3) + get(mol,'n3')/2+1,...
%        ng, ng);

igpadd = igcol; % line 142
qgp    = G_index(igpadd, :) + qk(qindex, :);

for igrow = 1:ng
  igadd = igrow; % line 150, for corresponding relations between G to q+G.
  Omega2(igadd) = 0;
  gg = G_index(igadd, :) - G_index(igpadd, :);
%  gg1D = mill2nl(gg, mol.n1, mol.n2, mol.n3);
  
%   disp([num2str(gg), ' ',num2str(gg1D),  ' ', num2str(gg1Dstar)])
  kadd = findvector(gg, gvec, setting_method);
%   kadd = find(idxnz == gg1D);
%   if isempty(kadd)
%     continue;
%   end
	if kadd == 0; continue; end
	rho_g_minus_gp = rho(kadd); %% sum over kpoints if necessary.
  if (igadd ~= 1 || q_is_not_zero)
  %  Omega2(igrow) =  (qgp * bdot * precalc(igrow, :)') * rho_g_minus_gp;     
    Omega2(igadd) =  (qgp * GWinfo.bdot * precalc(igrow, :)') * rho_g_minus_gp;
  elseif igpadd == 1
    Omega2(igadd) = fact_wpeff * rho_g_minus_gp;
  else
    Omega2(igadd) = 0;
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
%       Omega2(igrow) =  (qgp * GWinfo.bdot * precalc(igrow, :)') * rho_g_minus_gp;
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
%       Omega2(igrow) =  (qgp * GWinfo.bdot * precalc(igrow, :)') * rho_g_minus_gp;
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
