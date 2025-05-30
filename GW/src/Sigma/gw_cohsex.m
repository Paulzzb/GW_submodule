function [Esx_x, Ecoh] = gw_cohsex(GWinfo, options)
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
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
Dcoul(1,1) = GWinfo.coulG0;
gvec = GWinfo.gvec;
ev    = GWinfo.ev * ry2ev;
Dcoul = Dcoul * ry2ev;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start calculation
if (options.ISDFCauchy.isISDF)
  startISDF = tic;
  optISDF = options.ISDFCauchy;
  Phig = GWinfo.Z;
  nr = gvec.nfftgridpts; nb = size(Phig, 2); 
  psir = zeros(nr, nb);
  for iband = 1:nb
    fftbox1 = put_into_fftbox(Phig(:, iband), gvec.idxnz, gvec.fftgrid);
    fftbox1 = nr / vol * do_FFT(fftbox1, gvec.fftgrid, 1);
    psir(:, iband) = reshape(fftbox1, nr, []);
  end
  
  vsrank_mu = ceil(sqrt(nv_oper*n_ener)  * optISDF.vsrank_ratio);
  optISDF.isdfoptions.rank = vsrank_mu;
  [vsind_mu, vsgzeta_mu] = isdf_main('vs', Phig, nv-nv_oper+1:nv, ...
          nv-nv_ener+1:nv+nc_ener, gvec, vol, optISDF);
  vsgzeta_mu = conj(vsgzeta_mu);
  
  % In this case, ISDF for vc, vs, ss are all needed.
  vcrank_mu = ceil(sqrt(nv_oper*nc_oper) * optISDF.vcrank_ratio);
  ssrank_mu = ceil(sqrt(n_ener *n_ener)  * options.ISDFCauchy.ssrank_ratio);
  optISDF.isdfoptions.rank = vcrank_mu;
  [vcind_mu, vcgzeta_mu] = isdf_main('vc', Phig, nv-nv_oper+1:nv, ...
          nv+1:nv+nc_oper, gvec, vol, optISDF);
  vcgzeta_mu = conj(vcgzeta_mu);
  
  optISDF.isdfoptions.rank = ssrank_mu;
  [ssind_mu, ssgzeta_mu] = isdf_main('ss', Phig, nv-nv_oper+1:nv+nc_oper, ...
          nv-nv_ener+1:nv+nc_ener, gvec, vol, optISDF);
  ssgzeta_mu = conj(ssgzeta_mu);
  
  psirvc = psir(vcind_mu, :);
  psirvs = psir(vsind_mu, :);
  psirss = psir(ssind_mu, :);

  timeforISDF = toc(startISDF);



  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Manipulate operators

  % Calculate C*Omega^(-1)*C', use direct method or cauchy integral.
  startinveps = tic; 

  if options.ISDFCauchy.isCauchy
    Phi = psirvc(:, nv-nv_oper+1:nv); 
    Psi = psirvc(:, nv+1:nc_oper+nv); 
    evOcc = ev(nv-nv_oper+1:nv);
    evUnocc = ev(nv+1:nv+nc_oper);
    
    [COmegaCresult, ~, ~] = ...
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
  epsg_main = inv(COmegaCresult) - vcgzeta_mu' * Dcoul * vcgzeta_mu; 
  timeforinveps = toc(startinveps);
  fprintf('Time for inveps = %.4f.\n', timeforinveps);

  
  startSigma = tic; 
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
  timeForSigma = toc(startSigma);
  fprintf('Sigma operator time = %f.\n', timeForSigma);

  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  % Calculate self-energies

  timeForSelfE = tic;
  Psivs = conj(psirvs(:, nv-nv_ener+1:nv+nc_ener)); 
  Psiss = conj(psirss(:, nv-nv_ener+1:nv+nc_ener)); 
  
  Esx_x = - Psivs' * Sigma_sex_x * Psivs / GWinfo.vol;
  Ecoh   = 0.5 * Psiss' * Sigma_coh   * Psiss / GWinfo.vol;

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
  timeforW = toc(startforW);
  fprintf('timeforW = %.4f sec.\n', timeforW); 

  startEsx_xExEch = tic;
  Esx_x  = zeros(n_ener);
  Ecoh    = zeros(n_ener);

  for ioper = nv-nv_oper+1:nv
    Mgvn = mtxel_sigma(ioper, GWinfo, options.Groundstate, (nv-nv_ener+1:nv+nc_ener));
    Mgvn = conj(Mgvn);
    W1Mgvn = Dcoul * Mgvn;

    W1Mgvn = inveps \ W1Mgvn;
    W1Mgvn = W1Mgvn - Dcoul * Mgvn;
    Esx_x = Esx_x - Mgvn' * W1Mgvn / GWinfo.vol;
    Ecoh = Ecoh + 0.5 * Mgvn' * W1Mgvn / GWinfo.vol;
  end
  clear Mgvn W1Mgvn;

  for ioper = nv+1:nv+nc_oper
    Mgcn = mtxel_sigma(ioper, GWinfo, options.Groundstate, (nv-nv_ener+1:nv+nc_ener));
    Mgcn = conj(Mgcn);
    W1Mgcn = Dcoul * Mgcn;
    W1Mgcn = inveps \ W1Mgcn;
    W1Mgcn = W1Mgcn - Dcoul * Mgcn;
    
    Ecoh = Ecoh + 0.5 * Mgcn' * W1Mgcn / GWinfo.vol;
  end
  clear W1Mgcn Mgcn;

  timeforEsx_xExEch = toc(startEsx_xExEch);
  fprintf('Time for Ex, Esx_x and finite cutoff Ecoh = %.4f.\n', ...
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
end


Esx_x = real(diag(Esx_x));
Ecoh = real(diag(Ecoh));



% save(['GWenergy_',num2str(vcrank_mu),'_',num2str(vsrzeta_mu),'_',num2str(ssrank_mu),'.mat'], GWenergy);
% save(['GWenergy', name_of_mol, '_',num2str(vcrank_mu), '_', num2str(vsrank_mu), '_', num2str(ssrank_mu), '.mat'], 'GWenergy');
% save(['timeinfo', name_of_mol], 'time*');
end
