function gw_cohsex_gpu(GWinfo, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is used to calculate GW quasiparticle energies under COHSEX approximation.
%
% Input:
% GWinfo: A structure that contains ground state calculation results
% and parameters of the molecules.
% Check gwsetup.m or README for details.
%
% options: 
%   isISDF:    true or false.
%   iscauchy:  use Cauchy integral or not, only applicable if isISDF = True.
%   optionscauchy: options when doing cauchy integral, check introduction of 
%                   COmegaCstar.m for details.
%   vcrank_mu: ISDF coefficient for valence orbitals and conduction orbitals.  
%   vsrank_mu: ISDF coefficient for valence orbitals and all orbitals.  
%   ssrank_mu: ISDF coefficient for all orbitals and all orbitals.  
%   (*rank_mu is only applicable if isISDF = true.)
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialization
% global nv nv_ener nv_oper nc nc_ener nc_oper n_oper n_ener
% global ng nr ne vol bdot
% global ry2ev
% global im
% global TOL_SMALL
% ZERO = 0.0;
% CZERO = 0.0 + sqrt(-1) * 0.0;
startGW = tic;
nameConstants = fieldnames(options.Constant);
for i = 1:numel(nameConstants)
    fieldname = nameConstants{i};
    value = options.Constant.(fieldname);    
    if ~isempty(value)
      strEval = sprintf('%s = %.16f;', fieldname, value);
      eval(strEval);
    end
end

GWinfo.Z = GWinfo.Z * sqrt(GWinfo.vol); % Turn to \int_V \abs{psi(r)}^2 dr = 1.
Z     = GWinfo.Z;
ev    = GWinfo.ev;
Vxc   = GWinfo.Vxc;
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
gvec = GWinfo.gvec;
Dcoul(1,1) = GWinfo.coulG0;
fprintf('nr = %d, ng = %d, n_oper = %d\n', nr, ng, n_oper);



% Normalize the wavefunc in Fourier space.
startC2R = tic;
psir = zeros(nr, nv+nc);
for iband = 1:nv+nc
  fftbox1 = put_into_fftbox(Z(:, iband), gvec.idxnz, gvec.fftgrid);
  fftbox1 = gvec.nfftgridpts / GWinfo.vol * do_FFT(fftbox1, gvec.fftgrid, 1);
  psir(:, iband) = reshape(fftbox1, gvec.nfftgridpts, []);
end
timeforC2R = toc(startC2R);
fprintf('Time for psir = %.4f.\n', timeforC2R)
clear Z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start calculation
if (options.ISDFCauchy.isISDF)
  vcrank_mu = ceil(sqrt(nv_oper*nc_oper) * options.ISDFCauchy.vcrank_ratio);
  vsrank_mu = ceil(sqrt(nv_oper*n_ener)  * options.ISDFCauchy.vsrank_ratio);
  ssrank_mu = ceil(sqrt(n_ener *n_ener)  * options.ISDFCauchy.ssrank_ratio);
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Start ISDF
  % Generally, we have relations 
  %     Mr{xx} = xxrzeta_mu * Mrxx_mu;
  vcgzeta_mu = zeros(ng, vcrank_mu);
  optionsISDF = options.ISDFCauchy;


  startVC = tic;
  optionsISDF.isdfoptions.rank = vcrank_mu;
  vcind_mu = isdf_indices(conj(psir(:,nv-nv_oper+1:nv)), ...
                               (psir(:,nv+1:nv+nc_oper)), optionsISDF);
  vcgzeta_mu = isdf_kernelg(conj(psir(:,nv-nv_oper+1:nv)), ...
              (psir(:,nv+1:nv+nc_oper)), vcind_mu, gvec, vol);
  vcgzeta_mu = conj(vcgzeta_mu);
  timeforVC = toc(startVC);
  

  startVS = tic;
  vsgzeta_mu = zeros(ng, vsrank_mu);
  optionsISDF.isdfoptions.rank = vsrank_mu;
  vsind_mu = isdf_indices(conj(psir(:,nv-nv_oper+1:nv)), ...
                                (psir(:,nv-nv_ener+1:nv+nc_ener)), optionsISDF);
  vsgzeta_mu = isdf_kernelg(conj(psir(:,nv-nv_oper+1:nv)), ...
              (psir(:,nv-nv_ener+1:nv+nc_ener)), vsind_mu, gvec, vol);
  vsgzeta_mu = conj(vsgzeta_mu);
  timeforVS = toc(startVS);


  startSS = tic;
  ssgzeta_mu = zeros(ng, ssrank_mu);
  optionsISDF.isdfoptions.rank = ssrank_mu;
  ssind_mu = isdf_indices(conj(psir(:,nv-nv_oper+1:nv+nc_oper)), ...
                                (psir(:,nv-nv_ener+1:nv+nc_ener)), optionsISDF);
  ssgzeta_mu = isdf_kernelg(conj(psir(:,nv-nv_oper+1:nv+nc_oper)), ...
              (psir(:,nv-nv_ener+1:nv+nc_ener)), ssind_mu, gvec, vol);
  ssgzeta_mu = conj(ssgzeta_mu);
  timeforSS = toc(startSS);

  % startaqsch = tic;
  % aqsch = zeros(ng, n_ener);
  % tmpC = psir(ssind_mu, nv-nv_ener+1:nv+nc_ener);
  % tmpC = abs(tmpC).^2;
  % aqsch = ssgzeta_mu * tmpC;
  % clear ssrzeta_mu;
  % timeforaqsch = toc(startaqsch);
  % clear aqsch

  psirvc = psir(vcind_mu, :);
  psirvs = psir(vsind_mu, :);
  psirss = psir(ssind_mu, :);
  clear psir;

  fprintf('ISDF time = %.4f.\n', timeforVC+timeforVS+timeforSS);
  % fprintf('Aqsch time = %.4f.\n', timeforaqsch);
  % ISDF DONE

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Manipulate operators

  % Calculate C*Omega^(-1)*C', use direct method or cauchy integral.
  startinveps = tic; 

  if options.ISDFCauchy.isCauchy
    Phi = psirvc(:, nv-nv_oper+1:nv); 
    Psi = psirvc(:, nv+1:nc_oper+nv); 
    evOcc = ev(nv-nv_oper+1:nv);
    evUnocc = ev(nv+1:nv+nc_oper);
    
    [COmegaCresult, time_c, relError] = ...
          COmegaCstar(Phi, Psi, evOcc, evUnocc, options.ISDFCauchy.optionsCauchy);
    clear Phi Psi evOcc evUnocc;
    COmegaCresult = COmegaCresult * 4;
  else
    Mgvc = zeros(vcrank_mu, nc_oper);
    COmegaCresult = zeros(vcrank_mu, vcrank_mu);
    scal = 4.0;
  %  scal = 4.0 / GWinfo.vol; % for kpoints & spin, change it according to
                % epsilon_main.f90 and chi_summation.f90
    for ind_nv = nv-nv_oper+1:nv
      Mgvc = conj(psirvc(:, ind_nv)) .* psirvc(:, nv+1:nv+nc_oper);
      Mgvc = conj(Mgvc);
      eden = 1 ./ (ev(ind_nv) - ev(nv+1:nv+nc_oper));
      COmegaCresult = COmegaCresult + scal * Mgvc * diag(eden) * Mgvc';
    end
  end
  COmegaCresult = COmegaCresult / GWinfo.vol;
  % Calculate ( (C*Omega*C)^(-1) ...
  %                               - vcgzeta_mu' * Dcoul * vcgzeta_mu)^(-1) 
  epsg_main = inv(COmegaCresult) - vcgzeta_mu' * Dcoul * vcgzeta_mu; 
  timeforinveps = toc(startinveps);
  fprintf('Time for inveps = %.4f.\n', timeforinveps);

  
  startSigma = tic; 
  Dcoul(1,1) = GWinfo.coulG0;
  % Calculate operator \Sigma_{COH}, \Sigma_{SEX_X}, and \Sigma_{X} formally.
  Phivs = psirvs(:, nv-nv_oper+1:nv); 
  Phiss = psirss(:, nv-nv_oper+1:nv+nc_oper); 
  epsvc = vcgzeta_mu / epsg_main;
  epsvcDcoulvs = epsvc' * Dcoul * vsgzeta_mu;
  epsvcDcoulss = epsvc' * Dcoul * ssgzeta_mu;
  vcDcoulvs = vcgzeta_mu' * Dcoul * vsgzeta_mu;
  vcDcoulss = vcgzeta_mu' * Dcoul * ssgzeta_mu;
  W1_mu = vcDcoulvs' * (epsvcDcoulvs);
  W1_mu_1 = vcDcoulss' * (epsvcDcoulss);
  Sigma_sex_x = W1_mu .* (Phivs * Phivs'); 
  Sigma_coh   = W1_mu_1 .* (Phiss * Phiss');
  clear W1_mu vcDcoulvs vcDcoulss epsvcDcoulvs epsvcDcoulss;
  vsDcoulvs = vsgzeta_mu' * Dcoul * vsgzeta_mu;
  Sigma_x     = vsDcoulvs .* conj(Phivs * Phivs'); 
  clear vsDcoulvs;
  timeForSigma = toc(startSigma);
  fprintf('Sigma operator time = %f.\n', timeForSigma);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate self-energies

  timeForSelfE = tic;
  Psivs = conj(psirvs(:, nv-nv_ener+1:nv+nc_ener)); 
  Psiss = conj(psirss(:, nv-nv_ener+1:nv+nc_ener)); 
  
  Esx_x = Psivs' * Sigma_sex_x * Psivs / GWinfo.vol;
  Ech   = Psiss' * Sigma_coh   * Psiss / GWinfo.vol;
  Ex    = Psivs' * Sigma_x     * Psivs / GWinfo.vol;

  clear Psivs Psiss;
  timeforSelfE = toc(timeForSelfE);
  fprintf('Self-energies time = %f.\n', timeforSelfE);
  
%%%%%%%%%%%%%
  % This part is added for debugging purpose.
  % Ex_2 = zeros(n_ener, n_ener);
  % vsDcoulvs = vsgzeta_mu' * Dcoul * vsgzeta_mu;
  % for ioper = nv-nv_oper+1:nv
  %   Mgvn = conj(psir(vsind_mu, ioper)) .* psir(vsind_mu, nv-nv_ener+1:nv+nc_ener);
  %   Mgvn = (Mgvn);
  %   Ex_2 = Ex_2 + Mgvn' * vsDcoulvs * Mgvn / GWinfo.vol;
  % end
  % clear Mgcn;

%%%%%%%%%%%%

%   % Calculate exactCH
%   startExactCH = tic;
%   L = W1_mu_1;
%   % L = chol(W1_mu_1);
%   tmpC = zeros(ssrank_mu, n_ener);
%   aqsch = zeros(nr, n_ener);
%   Wrr = zeros(nr, 1); 
%   tmpssrzetamu = zeros(nr, 1);
% 
%   tmpC = psir(ssind_mu, nv-nv_ener+1:nv+nc_ener);
%   tmpC = abs(tmpC).^2;
%   ssrzetamu_ = zeros(nr, ssrank_mu);
%   ssrzetamu_right = zeros(nr, ssrank_mu);
% 
%   for indss = 1:ssrank_mu
%     fftbox1 = put_into_fftbox(ssgzeta_mu(:, indss), gvec.idxnz, gvec.fftgrid);
%     fftbox1 = do_FFT(fftbox1, gvec.fftgrid, 1) * GWinfo.vol;
%     tmpssrzetamu = reshape(fftbox1, gvec.nfftgridpts, []);
%     aqsch = aqsch + tmpssrzetamu * tmpC(indss, :);
%     ssrzetamu_ = ssrzetamu_ + tmpssrzetamu * L(indss, :);
%     ssrzetamu_right(:, indss) = tmpssrzetamu;
%   end
% 
%   for ir = 1:nr
%     Wrr(ir) = ssrzetamu_(ir, :) * ssrzetamu_right(ir, :)';
%   end
%     
%   ach = sum(aqsch .* Wrr, 2) / 2;
%   
%   ach * ry2ev;
% 
%   timeForExactCH = toc(startExactCH)
  
  % Calculate exactCH
  % startExactCH = tic;  
  % achx = zeros(n_ener, 1);

  % % We already have aqsch. 
  % W1 = Dcoul * vcgzeta_mu * invepsg_main * vcgzeta_mu' * Dcoul;
  % for igcol = 1:ng
  %   gpp_list = gvec.components - gvec.components(igcol, :);
  %   [iout_list, gindex_list] = findvector(gpp_list, gvec);
  %   achx = achx + aqsch(iout_list(gindex_list), :).' * W1(gindex_list, igcol) / 2 / GWinfo.vol;
  % end % for igcol = 1:ng
  %  timeForExactCH = toc(startExactCH);
  % fprintf('Time for Exact CH = %f.\n', timeForExactCH)


  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate energies
  GWenergy.ev = ev(nv-nv_ener+1:nv+nc_ener) * ry2ev;
  GWenergy.Ex = real(diag(0.0 - Ex) * ry2ev);
  GWenergy.Esx_x = real(diag(0.0 - Esx_x) * ry2ev);
  GWenergy.Ech = real(diag(0.5 * Ech * ry2ev));
  % GWenergy.Ech_exact = real(achx * ry2ev);
  GWenergy.Sig = real(GWenergy.Ex + GWenergy.Esx_x + GWenergy.Ech);
  GWenergy.Vxc = Vxc(nv-nv_ener+1:nv+nc_ener) * ry2ev;
  GWenergy.eqp = GWenergy.ev - GWenergy.Vxc + GWenergy.Sig;

else % No ISDF version
%   gvec_bgw = readmatrix('gvec_components.csv');
%   ks2bgw = index_ks2bgw(GWinfo.coulG(:, 1:3), gvec_bgw(1:81, :));
%   [~, bgw2ks] = sort(ks2bgw);
 
  startforW = tic;

  Mgvc = zeros(ng, nc_oper);
  inveps = zeros(ng, ng);
  scal = 4.0;
%  scal = 4.0 / GWinfo.vol; % for kpoints & spin, change it according to
              % epsilon_main.f90 and chi_summation.f90
  for ind_nv = nv-nv_oper+1:nv
    Mgvc = mtxel_sigma(ind_nv, GWinfo, options.Groundstate, (nv+1:nv+nc_oper));
    Mgvc = conj(Mgvc);
    eden = 1 ./ (ev(ind_nv) - ev(nv+1:nv+nc_oper));
    inveps = inveps + scal * Mgvc * diag(eden) * Mgvc';
  end
  inveps = eye(ng) - Dcoul * inveps / vol;
  % Real inveps
  % inveps = inv(inveps);
%  inveps = inveps - eye(ng);  



  % Sigma part
  Dcoul(1) = GWinfo.coulG0;
  % W1     = inveps * Dcoul;
  % W1     = W1 - Dcoul;
  % W1_     = (inveps + eye(ng) * Dcoul;
  timeforW = toc(startforW);
  fprintf('timeforW = %.4f sec.\n', timeforW); 

  startEsx_xExEch = tic;
  Esx_x  = zeros(n_ener);
  Ex     = zeros(n_ener);
  Ech    = zeros(n_ener);

  for ioper = nv-nv_oper+1:nv
    Mgvn = mtxel_sigma(ioper, GWinfo, options.Groundstate, (nv-nv_ener+1:nv+nc_ener));
    Mgvn = conj(Mgvn);
    W1Mgvn = Dcoul * Mgvn;
    Ex = Ex + Mgvn' * W1Mgvn / GWinfo.vol;

    W1Mgvn = inveps \ W1Mgvn;
    W1Mgvn = W1Mgvn - Dcoul * Mgvn;
    Esx_x = Esx_x + Mgvn' * W1Mgvn / GWinfo.vol;
    Ech = Ech + Mgvn' * W1Mgvn / GWinfo.vol;
  end
  clear Mgvn W1Mgvn;

  for ioper = nv+1:nv+nc_oper
    Mgcn = mtxel_sigma(ioper, GWinfo, options.Groundstate, (nv-nv_ener+1:nv+nc_ener));
    Mgcn = conj(Mgcn);
    W1Mgcn = Dcoul * Mgcn;
    W1Mgcn = inveps \ W1Mgcn;
    W1Mgcn = W1Mgcn - Dcoul * Mgcn;
    
    Ech = Ech + Mgcn' * W1Mgcn / GWinfo.vol;
  end
  clear W1Mgcn Mgcn;

  timeforEsx_xExEch = toc(startEsx_xExEch);
  fprintf('Time for Ex, Esx_x and finite cutoff Ech = %.4f.\n', ...
          timeforEsx_xExEch);  

  % timeforEchExact = tic;
  % gvec = GWinfo.gvec;
  % isrtrq = 1:gvec.ng;
  % isrtrqi = 1:gvec.ng;
  % ncoulch = ng;
  % achx = zeros(n_ener, 1);
  
  
  % Another implememtation of exact_ch
  % 1. transform W(g1, g2) --> W(r, g2);
  % 2. calculate W(r, r) elementwise.
  % 3. calculate <nn|diag(W)>.
  % startexactch1 = tic;
%  startaqsch = tic;
%  aqsch = zeros(nr, n_ener);
%  fftbox1 = zeros(gvec.fftgrid);
%  for iener = nv-nv_ener+1:nv+nc_ener
%    fftbox1 = put_into_fftbox(Z(:, iener), gvec.idxnz, gvec.fftgrid);
%    fftbox1 = gvec.nfftgridpts / GWinfo.vol * do_FFT(fftbox1, gvec.fftgrid, 1);
%    fftbox1 = fftbox1 .* conj(fftbox1);
%    switch lower(options.Groundstate.setting_method)
%      case 'kssolv'
%        aqsch(:, iener) =  reshape(fftbox1, [], 1); 
%%      case 'bgw'
%%        
%%        aqsch()
%      otherwise
%        fprintf('not support!\n')
%        error()
%    end
%  end
%  timeaqsch = toc(startaqsch)
%  startWrr = tic;
%  W_rg_tmp = zeros(nr, 1);
%  tmparray = zeros(ng, 1);
%  tmparrayr = zeros(nr, 1);
%  % Construct rindex to generate DFT matrix explicitly.
%  mol = GWinfo.mol;
%  n1 = mol.n1;
%  n2 = mol.n2;
%  n3 = mol.n3;
%  rindex = ones(3, n1, n2, n3);
%  for i = 1:n1
%    rindex(1, i, :, :) = i-1;
%  end
%  for i = 1:n2
%    rindex(2, :, i, :) = i-1;
%  end
%  for i = 1:n3
%    rindex(3, :, :, i) = i-1;
%  end
%
%  rindex = reshape(rindex, 3, nr);
%  rindex = rindex.';
%
%  for igcol = 1 : ng
%    fftbox1 = put_into_fftbox(W1(:, igcol) * GWinfo.vol, gvec.idxnz, gvec.fftgrid);
%    fftbox1 = gvec.nfftgridpts / GWinfo.vol * do_FFT(fftbox1, gvec.fftgrid, 1);
%    W_rg_tmp = reshape(fftbox1, [], 1);
%%     if (igcol ~= 1)
%%       tmparray(igcol-1) = 0;
%%     end
%%     tmparray(igcol) = 1;
%%     fftbox1 = put_into_fftbox(tmparray, gvec.idxnz, gvec.fftgrid);
%%     fftbox1 = do_FFT(fftbox1, gvec.fftgrid, -1);
%%     tmparrayr = reshape(fftbox1, gvec.nfftgridpts, []);
%    
%    gindex = gvec.components(igcol, :) ./ ...
%             [GWinfo.mol.n1, GWinfo.mol.n2, GWinfo.mol.n3] * 2 * pi;
%    tmparrayr_ = exp(-sqrt(-1) * rindex * gindex');
%%    if norm(tmparrayr_ - tmparrayr) >= 1e-10
%%      fprintf('tmparrayr ~= tmparrayr_.\n');
%%      error()
%%    end
%
%    W_rg_tmp = W_rg_tmp .*  tmparrayr_ / vol;
%    achx = achx + (sum(W_rg_tmp .* aqsch, 1))' * GWinfo.vol / gvec.nfftgridpts;
%  end
%  timeWrr = toc(startWrr)
  % timeexactch1 = toc(startexactch1);
  % fprintf('Time of exactch1 = %f.\n', timeexactch1);
  % achx * ry2ev / 2;
  % achx_ = - achx / 2;
  % startexactch2 = tic;
  % achx = zeros(n_ener, 1);
  % W1 = inveps - eye(ng); 
  % for iener = nv-nv_ener+1:nv+nc_ener
%  for iener = 1:1
%    if aqsFlag
%      aqsch = GWinfo.aqs{iener}(:, iener);
%    else
%      aqsch = mtxel_sigma(iener, GWinfo, options.Groundstate, iener);
%    end

%     for igcol = 1 : 20
%       achxtemp_gp = CZERO;
%       gpp_list = gvec.components - gvec.components(igcol, :);
%       iout2_list = findvector(gpp_list, gvec);
%       for igrow = 1 : ng
%         I_epsggp = W1(igrow, igcol);
%          
% %         if abs(I_epsggp) < TOL_SMALL
% %           continue; 
% %         end
%       
%         if igrow ~= igcol
% %           gpp = gvec.components(isrtrq(igrow), :) - gvec.components(isrtrq(igcol), :);
% %           iout = findvector(gpp, gvec);
%       
%           iout = iout_list(igrow);
%           if iout ~= iout2
%             error()
%           end
%           if iout == 0
%             continue; % Equivalent to Fortran's "cycle"
%           end
%       
%           igpp = isrtrqi(iout);
%       
%           if igpp < 1 || igpp > ncoulch
%             continue; % Equivalent to Fortran's "cycle"
%           end
%       
% %           gpp = gvec.components(isrtrq(igcol), :) - gvec.components(isrtrq(igrow), :);
% %           iout = findvector(gpp, gvec);
% %           % iout2 = findvector_swp(gpp, gvec);
% %       
% %           if iout == 0; continue; end
% %       
% %           igpp2 = isrtrqi(iout);
% %       
% %           if (igpp2 < 1 || igpp2 > ncoulch); continue; end
%         else
%           iout = iout_list(igrow);
% %           gpp = zeros(1, 3);
% %           iout = findvector(gpp, gvec);
% %           iout2 = findvector_swp(gpp, gvec);
% %           if iout ~= iout2
% %             error()
% %           end
%       
%           if iout == 0
%             continue; % Equivalent to Fortran's "cycle"
%           end
%       
%           igpp = isrtrqi(iout);
%       
%           if igpp < 1 || igpp > ncoulch
%             continue; %
%           end
%         end
%       
%         schx = aqsch(igpp) * I_epsggp;
%         achxtemp_gp = achxtemp_gp + schx;
%       end % for igrow
%       achx(iener) = achx(iener) + achxtemp_gp * vcoul(igcol) * 0.5 / GWinfo.vol;
%     end % for igcol
%   end % for iener
% Exact CH   
  
  % achx = zeros(n_ener, 1);
  % aqsch = zeros(ng, n_ener);
  % for iener = nv-nv_ener+1 : nv+nc_ener
  %   if aqsFlag
  %     aqsch(:, iener-nv+nv_ener) = GWinfo.aqs{iener}(:, iener);
  %   else
  %     aqsch(:, iener-nv+nv_ener) = mtxel_sigma(iener, GWinfo, options.Groundstate, iener);
  %   end
  % end
  
  % W1 = W1 * Dcoul;
  % for igcol = 1:ng
  %   gpp_list = gvec.components - gvec.components(igcol, :);
  %   [iout_list, gindex_list] = findvector(gpp_list, gvec);
  %   achx = achx + aqsch(iout_list(gindex_list), :).' * W1(gindex_list, igcol) / 2 / GWinfo.vol;
  % end % for igcol = 1:ng
  
  

  % timeexactch2 = toc(startexactch2);
  % fprintf('Time for exactch 2 = %f.\n', timeexactch2);
  % achx * ry2ev;
  % Calculate energies
  GWenergy.ev = ev(nv-nv_ener+1:nv+nc_ener) * ry2ev;
  GWenergy.Ex = real(diag(0.0 - Ex) * ry2ev);
  GWenergy.Esx_x = real(diag(0.0 - Esx_x) * ry2ev);
  GWenergy.Ech = real((diag(0.5 * Ech))* ry2ev);
  % GWenergy.Ech = real(achx) * ry2ev;
  % GWenergy.Ech = real(achx)* ry2ev;
% GWenergy.Ech1 = real((diag(0.5 * Ech) + achx)* ry2ev);
%  GWenergy.Ech
  GWenergy.Sig = real(GWenergy.Ex + GWenergy.Esx_x + GWenergy.Ech);
  GWenergy.Vxc = Vxc(nv-nv_ener+1:nv+nc_ener) * ry2ev;
  GWenergy.eqp = GWenergy.ev - GWenergy.Vxc + GWenergy.Sig;

end





% save(['GWenergy_',num2str(vcrank_mu),'_',num2str(vsrzeta_mu),'_',num2str(ssrank_mu),'.mat'], GWenergy);
% save(['GWenergy', name_of_mol, '_',num2str(vcrank_mu), '_', num2str(vsrank_mu), '_', num2str(ssrank_mu), '.mat'], 'GWenergy');
% save(['timeinfo', name_of_mol], 'time*');
timeforGW = toc(startGW);
save(options.GWCal.fileName, 'GWenergy', 'time*') 
fprintf('Time for GW Calculation under COHSEX approximation = %f.\n', timeforGW);
fprintf('Result is saved in %s.\n', options.GWCal.fileName)
end
