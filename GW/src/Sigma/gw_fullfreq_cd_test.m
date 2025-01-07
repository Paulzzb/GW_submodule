function [Eqp0,Eqp1]=gw_fullfreq_cd(GWinfo, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% GWinfo is a structure that contains ground state calculation results
% and parameters of the molecules
%   nv      --- number of valence states
%   Z       --- contains the eigenvecgtors from the KSDFT calculation or Quantum-Espresso calculation
%   ev      --- contains the corresponding eigenvalues
%   vol     --- volume of the unit cell
%   ntot    --- total number of grid points on which the wavefunction is sampled
% 
% eta     --- Lorentzian broadening factor that turns a Dirac delta into a 
%             a smoother peak
% nv_oper --- the number of valence bands (KS orbitals) from the Fermi 
%             level included in the kernel calculation. 
% nc_oper --- the number of empty bands (KS orbitals) included in 
%             the kernel calculation
% n_oper  --- the number of empty bands (KS orbitals) included in
%             the Coulomb-hole (Ech) calculation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


testInput = true;
testInput = false;
testflag1 = true;
% testflag1 = false;
testflag2 = true;
testflag2 = false;
testflag3 = true;
testflag3 = false;


% Initialization
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

if testInput
  GWinfor2 = load("SiGWinfo.mat");
  GWinfor2 = GWinfor2.GWinfo;
  GWinfo.Z = GWinfor2.Z;
  GWinfo.ev = GWinfor2.ev;
  % GWinfo.aqs = GWinfor2.aqs;
end
% Initialization
GWinfo.Z     = GWinfo.Z * sqrt(vol);
Z     = GWinfo.Z;
ev    = GWinfo.ev;
F     = GWinfo.F;
Vxc   = GWinfo.Vxc;
nfreq_real = options.GWCal.nFreq - options.GWCal.nfreq_imag;
nfreq_imag = options.GWCal.nfreq_imag;
% real_freq = GWinfo.real_freq;
% imag_freq = GWinfo.imag_freq;
% nfreqreal = size(real_freq,2);
% nfreqimag = size(imag_freq,2);
% nfreqeval = 2*floor((options.max_freq_eval+TOL_SMALL)/options.delta_freq_step)+1;
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
aqsFlag = ~isempty(GWinfo.aqs);

if isfield(options.GWCal, 'nFreq')
  nFreq = options.GWCal.nFreq;
else
  error('options does not have field nFreq. Please use gwCalculation.m to do the jobs.');
end

if isfield(options.GWCal, 'dFreqGrid')
  dFreqGrid = options.GWCal.dFreqGrid;
else
  error('options does not have field dFreqGrid. Please use gwCalculation.m to do the jobs.');
end

if isfield(options.GWCal, 'dFreqBrd')
  dFreqBrd = options.GWCal.dFreqBrd;
else
  error('options does not have field dFreqBrd. Please use gwCalculation.m to do the jobs.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start calculation

startGW = tic;

%------------------------------------------------------------------------------
% Calculate Epsilon

if 1
  epsilon = zeros(ng, ng, nFreq);
  for ind_nv = nv-nv_oper+1:nv
    if aqsFlag
      Mgvc = GWinfo.aqs{ind_nv}(:, nv+1:nv+nc_oper);
    else
      Mgvc = mtxel_sigma(ind_nv, GWinfo, options.Groundstate, ...
            (nv+1:nv+nc_oper));
    end
  	Mgvc = conj(Mgvc);
    Eden = ev(ind_nv) - ev(nv+1:nv+nc_oper);
    for ifreq = 1 : nFreq
  		edenDRtmp = 0.5 .* (1.0 ./ (Eden - (dFreqGrid(ifreq)+dFreqBrd(ifreq))) ...
  		+ 1.0 ./ (Eden + (dFreqGrid(ifreq)+dFreqBrd(ifreq))) );
  		epsilon(:, :, ifreq) = epsilon(:, :, ifreq) + 4 * Mgvc * (edenDRtmp .* Mgvc');
  	end % for ifreq
  end % for ind_nv
  
  for ifreq = 1:nFreq
  	epsilon(:, :, ifreq) = eye(ng) - ...
  	GWinfo.coulG(:, 4) .* epsilon(:, :, ifreq) / vol;
  end
  
  if testflag1
    fprintf("||epsilon||_F = %.3e.\n", norm(epsilon(:)));
    dlmwrite('epsilon_.csv', epsilon, 'precision', 16);
  end
  
  for ifreq = 1:nFreq
  	epsilon(:, :, ifreq) = inv(epsilon(:, :, ifreq));
  end
  
  if testflag1
    fprintf("||epsilon^-1||_F = %.3e.\n", norm(epsilon(:)));
  end
  
  for ifreq = 1:nFreq
  	epsilon(:, :, ifreq) = eye(ng) - epsilon(:, :, ifreq);
  end
  
end      


if 0
  Mgvc = zeros(ng, nv_oper * nc_oper);
  if ~isempty(GWinfo.aqs)
    for ind_nv = nv-nv_oper+1:nv
      ind_nv_true = ind_nv + (nv - nv_oper);
      Mgvc(:, (1:nc_oper) + (ind_nv_true-1) * (nc_oper)) ...
       = GWinfo.aqs{ind_nv}(:, nv+1:nv+nc_oper);
    end
  else
    aqs = cell(nv + nc, 1);
    aqstmp = zeros(ng, nv + nc);
  %  gvec = GWinfo.gvec;
    for ind_aqs = 1:nv+nc
      aqstmp = mtxel_sigma(ind_aqs, GWinfo, options.Groundstate, ...
          1:nv+nc);
      aqs{ind_aqs} = aqstmp;
    end
    GWinfo.aqs = aqs;
    for ind_nv = nv-nv_oper+1:nv
      ind_nv_true = ind_nv + (nv - nv_oper);
      Mgvc(:, (1:nc_oper) + (ind_nv_true-1) * (nc_oper)) ...
       = aqs{ind_nv}(:, nv+1:nv+nc_oper);
    end
    clear aqs scal
  end  
  Mgvc = conj(Mgvc);
  
  if testflag1
    fprintf("||Mgvc||_F = %.3e.\n", norm(Mgvc(:)));
  end
  
  % Calculate epsilon
  startEpsilon = tic;
  
  epsilon = zeros(ng, ng, nFreq);
  Eden = (kron(ev(nv-nv_oper+1:nv), ones(nc_oper,1)) ... 
                  - kron(ones(nv_oper, 1), ev(nv+1:nv+nc_oper)));
  
  for ifreq = 1 : nFreq 
    edenDRtmp = 0.5 .* ...
    ( 1.0 ./ (Eden - (dFreqGrid(ifreq) + dFreqBrd(ifreq))) ... 
    + 1.0 ./ (Eden + (dFreqGrid(ifreq) + dFreqBrd(ifreq))) );
    epsilon(:,:,ifreq) = eye(ng) - 4 * GWinfo.coulG(:,4) .* Mgvc * (edenDRtmp .* Mgvc') / vol;
  end
  clear Mgvc edenDRtmp;
  
  if testflag1
    fprintf("||epsilon||_F = %.3e.\n", norm(epsilon(:)));
    dlmwrite('epsilon.csv', epsilon, 'precision', 16);
  end
  
  timeforEpsilon = toc(startEpsilon)
  fprintf('\n');
  fprintf('Time for epsilon = %f.\n', timeforEpsilon);
  fprintf('\n');
  
  % Compare the results with outputs from BGW
  if exist('epsRDyn', 'file')
    epsRDyn_bgw = readfile('epsRDyn', '(%f, %f)');
    epsRDyn_bgw = complex(epsRDyn_bgw{1}, epsRDyn_bgw{2});
    epsRDyn_bgw = reshape(epsRDyn_bgw, ng, ng, nFreq);
    epsRDyn_bgw = epsRDyn_bgw(bgw2ks, bgw2ks, :);
  end
  
  % fprintf('\n');
  % for ifreq = 1 : nFreq
  %   fprintf('Frequency index: %d\n', ifreq);
  %   norm_diff = norm(epsRDyn_bgw(:, :, ifreq) - epsilon(:, :, ifreq), 'fro');
  %   fprintf('Norm difference = %f.\n', norm_diff);
  %   fprintf('\n')
  % end
  
  %------------------------------------------------------------------------------
  % Calculate W and install in epsilon
  
  startInv = tic;
  
  for ifreq = 1 : nFreq
    epsilon(:,:,ifreq) = inv(epsilon(:,:,ifreq));
  end
  
  timeforInv = toc(startInv)
  fprintf('\n');
  fprintf('Time for inverting epsilon = %f.\n', timeforInv);
  fprintf('\n');
  
  if testflag1
    fprintf("||inveps||_F = %.3e.\n", norm(epsilon(:)));
  end
  
  
  for ifreq = 1 : nFreq
     epsilon(:,:,ifreq) = eye(ng) - epsilon(:,:,ifreq);
  end 
end

I_eps_array = epsilon;
options.GWCal.dFreqBrd = options.GWCal.dFreqBrd * ry2ev;
options.GWCal.dFreqGrid = options.GWCal.dFreqGrid * ry2ev;
clear epsilon;


% -----------------------------------------------------------------------------
Dcoul(1, 1) = GWinfo.coulG0;
% Calculate <n|Sigma|n>
% We first calculate part contributed from residues (named as sres), 
% then calculate part contributed from integral (named as sint).
startSigma = tic;


Esx_x = CZERO * zeros(n_ener, options.GWCal.nfreqeval);
Ech = CZERO * zeros(n_ener, options.GWCal.nfreqeval);
Ex = CZERO * zeros(n_ener, 1);



% Calculate residues' contribution, all in real frequency
if options.GWCal.freq_grid_shift < 2
  freq0_array = options.GWCal.freqevalmin;
else 
% Doubt
  freq0_array = ev(nv-nv_ener+1:nv+nc_ener) * ry2ev ...
        - options.GWCal.freqevalstep * (options.GWCal.nfreqeval - 1) / 2;
end
if (options.GWCal.freq_dep == 2)
  freq0_array = freq0_array + TOL_SMALL;
end

iwlda_array = zeros(n_ener, 1);

for ibandouter = nv-nv_ener+1 : nv+nc_ener
  ibandouter_aqsntemp = ibandouter - nv + nv_ener;

  e_lk = ev(ibandouter) * ry2ev;
  freq0 = freq0_array(ibandouter_aqsntemp);
  if aqsFlag
    aqsntemp = GWinfo.aqs{ibandouter_aqsntemp};
  else
    aqsntemp = mtxel_sigma(ibandouter, GWinfo, options.Groundstate);
  end
  % aqsntemp = ...
  % GWinfo.aqs{ibandouter_aqsntemp};
  aqsntemp = conj(aqsntemp);
  asxt = CZERO;
  acht = CZERO;
  asxtemp = CZERO;
  achtemp = CZERO; 
  achstemp = CZERO;
  asxDtemp = CZERO * zeros(options.GWCal.nfreqeval, 1);
  achDtemp = CZERO * zeros(options.GWCal.nfreqeval, 1);

  for ibandinner = nv-nv_oper+1 : nv+nc_oper % line 1214 of mtxel_cor
     fprintf('ibandouter = %d, ibandinner = %d.\n', ibandouter, ibandinner);
%   pause(1);
    ibandinner_aqsntemp = ibandinner + nv_oper - nv;
    e_n1kq = ev(ibandinner_aqsntemp) * ry2ev;
    STATEMENT = (ibandinner <= nv); % refers to line 1221 
    if (STATEMENT) % we need adjustment of flagocc here!!!
      flagocc = true; occ = 1;
    else
      flagocc = false; occ = 0;
    end 
    Ex(ibandouter_aqsntemp) = Ex(ibandouter_aqsntemp) ...
                   + (aqsntemp(:, ibandinner_aqsntemp)' ...
                   * Dcoul * aqsntemp(:, ibandinner_aqsntemp)) ...
                   / GWinfo.vol * occ;
    
    
    achstemp = CZERO;
    sres = zeros(options.GWCal.nfreqeval, 1) * CZERO;
    sint = zeros(options.GWCal.nfreqeval, 1) * CZERO;
    % Start sigma_cd()
    diffmin = Inf;
    iwlda = -1;
    iw_array = ((0 : options.GWCal.nfreqeval-1) ...
                * options.GWCal.freqevalstep)';
    diff_array = abs(freq0 + iw_array - e_lk);
    [diffmin, iwlda] = min(diff_array);
    iwlda_array(ibandouter_aqsntemp) = iwlda;
    % for iw = 1:options.nfreqeval
    %   diff = abs(freq0 + (iw-1) * options.freqevalstep);
    %   if (diff < diffmin)
    %     diffmin = idff;
    %     iwlda = iw;
    %   end
    % end  
    % if iwlda < 0
    %   error('iwlda is not assigned a value!');
    % end

    wxi = freq0 - e_n1kq + iw_array;

    
    for iw = 1:options.GWCal.nfreqeval
      wx = wxi(iw);
      if ( (wx >= 0.0)  == flagocc)
        continue;
      end
      occ_sign = wx / abs(wx); 
      wx = abs(wx);
      ifreq = 0;

      for ijk = 1 : nfreq_real-1
        if (wx >= options.GWCal.dFreqGrid(ijk) && wx<options.GWCal.dFreqGrid(ijk+1))
          ifreq = ijk;
        end
      end

      if (ifreq == 0)
        ifreq = options.GWCal.nFreq + 3;
      end

      if (ifreq >= nfreq_real)
        fprintf('ibandouter = %d, ibandinner = %d, wx = %f, E_max = %f', ...
                ibandouter, ibandinner, wx, options.GWCal.dFreqGrid(nfreq_real));
        warning('Out of range occur!!')
        ifreq = nfreq_real - 1;
      end
      
      sres_omega = CZERO;
      if (nfreq_real > 1)
        fact1 = (options.GWCal.dFreqGrid(ifreq+1) - wx) ...
        / (options.GWCal.dFreqGrid(ifreq + 1) - options.GWCal.dFreqGrid(ifreq));
        fact2 = (-options.GWCal.dFreqGrid(ifreq) + wx) ...
        / (options.GWCal.dFreqGrid(ifreq + 1) - options.GWCal.dFreqGrid(ifreq));
      end
%      fprintf('fact1 = %f, fact2 = %f\n', fact1, fact2);

      if nfreq_real > 1
        sres_omega = sres_omega + ...
        aqsntemp(:, ibandinner_aqsntemp)'  ...
        * (fact1 * I_eps_array(:, :, ifreq) + fact2 * I_eps_array(:, :, ifreq+1)) ...
        * Dcoul * aqsntemp(:, ibandinner_aqsntemp) / GWinfo.vol;
      else
        sres_omega = sres_omega + aqsntemp(:, ibandinner_aqsntemp)' ...
                   * I_eps_array(:, :, 1) * Dcoul / GWinfo.vol * aqsntemp(:, ibandinner_aqsntemp);
      end 
      sres(iw) = sres(iw) - occ_sign * sres_omega;
%      for igcol = 1:ng
%        indigp = igp;
%        igp = igp;
%        sres_omega_acc = CZERO;
%        if (nfreq_real > 1)
%          for igrow = 1:ng
%          end
%        end
%
%      end %for igcol = 1:ng       
    end % for iw = 1:options.nfreqeval

    sW_imag_freqs = CZERO * zeros(options.GWCal.nfreq_imag, 1);
    for iw = 1:nfreq_imag
      sW_imag_freqs(iw) = sW_imag_freqs(iw) + aqsntemp(:, ibandinner_aqsntemp)' ...
                   * I_eps_array(:, :, iw+nfreq_real) * Dcoul / GWinfo.vol...
                   * aqsntemp(:, ibandinner_aqsntemp);
    end
    imag_freqs = CZERO * zeros(nfreq_imag, 1);
    imag_freqs = imag(options.GWCal.dFreqGrid(nfreq_real+1 : end) ...
                      + options.GWCal.dFreqBrd(nfreq_real + 1: end));

    switch options.GWCal.cd_int_method
      case 0
        for ifreq = 1 : options.GWCal.nfreq_imag 
          if (ifreq == 1)
            freqStart = imag_freqs(ifreq);
          else
            freqStart = (imag_freqs(ifreq - 1) + imag_freqs(ifreq)) * 0.5;
          end
          if (ifreq == options.GWCal.nfreq_imag)
            freqEnd = imag_freqs(end);
          else
            freqEnd = (imag_freqs(ifreq) + imag_freqs(ifreq + 1)) * 0.5;
          end
          sint = sint + sW_imag_freqs(ifreq) ...
                      * (atan(freqEnd ./ wxi) - atan(freqStart ./ wxi) );           
        end % for ifreq = 1: options.nfreq_imag
      % case 2
      % case 3
      otherwise
        error(['Options.cd_int_method = ', num2str(options.GWCal.cd_int_method), ' is not supported yet!']);
    end % switch options.cd_int_method

    asxDtemp = asxDtemp + sres;
    achDtemp = achDtemp + sint / pi;

  end %for ibandinner
  asxt = asxt + asxDtemp;
  acht = acht + achDtemp;
  Esx_x(ibandouter_aqsntemp, :) = asxt;
  Ech(ibandouter_aqsntemp, :)   = acht;
end% for ibandouter

Eres = Esx_x(:, (options.GWCal.nfreqeval + 1)/2)
Eint = Ech(:, (options.GWCal.nfreqeval + 1)/2)
Ex = real(Ex);

GWenergy = struct();
GWenergy.ev = ev(nv-nv_ener+1:nv+nc_ener) * ry2ev;
GWenergy.Ex = - Ex * ry2ev;
GWenergy.Eres = Eres * ry2ev;
GWenergy.Eint = Eint * ry2ev;
GWenergy.Sig = (GWenergy.Ex + GWenergy.Eres + GWenergy.Eint);
GWenergy.Vxc = Vxc(nv-nv_ener+1:nv+nc_ener) * ry2ev;
GWenergy.eqp = GWenergy.ev - GWenergy.Vxc + GWenergy.Sig;
save(options.GWCal.fileName, 'GWenergy') 

% iw = 1:nfreqeval;
% iw = repmat(iw,1,n_oper);
% ev_tmp = reshape((repmat(ev(1:nv+nc_oper)*ry2ev,1,nfreqeval))', ...
%                   1, n_oper*nfreqeval);
% psiphir = prod_states_gw(psir(:,nv-nv_ener+1:nv+nc_ener), ...
%                          conj(psir(:,nv-nv_oper+1:nv+nc_oper)));
% psiphi = F*psiphir;
% clear psiphir
% 
% wx = zeros(n_ener,n_oper*nfreqeval);
% occ_flag = zeros(n_ener,n_oper*nfreqeval);
% occ_sign = zeros(n_ener,n_oper*nfreqeval);
% fact1 = zeros(n_ener,n_oper*nfreqeval);
% fact2 = zeros(n_ener,n_oper*nfreqeval);
% ifreq = zeros(n_ener,n_oper*nfreqeval);
% 
% for i = 1:n_ener
%    
%    wx(i,:) = freq0(i) - ev_tmp + (iw-1) * options.delta_freq_step;
%    occ_flag(i,1:nv*nfreqeval) = 1;
%    tmp = find(xor(wx(i,:)>=0, occ_flag(i,:)==1));
%    %tmp_1 = tmp-(i-1)*n_oper*nfreqeval
%    occ_sign(i,tmp) = sign(wx(i,tmp));
%    wx_1 = zeros(1,n_oper*nfreqeval);
%    wx_1(tmp) = abs(wx(i,tmp));
%    for j = 1:nfreqreal-1
%       if(wx_1 >= real_freq(j)*ry2ev & wx_1 < real_freq(j+1)*ry2ev)
%          ifreq(i,tmp) = j;
%       end
%    end
%    fact1(i,tmp) = (real_freq(ifreq(i,tmp)+1) * ry2ev - wx_1(tmp)) ./ ...
%       (real_freq(ifreq(i,tmp)+1) - real_freq(ifreq(i,tmp))) / ry2ev;
%    fact2(i,tmp) = (wx_1(tmp) - real_freq(ifreq(i,tmp)) ) ./ ...
%       (real_freq(ifreq(i,tmp)+1) - real_freq(ifreq(i,tmp))) / ry2ev;
% end
% 
% sres = zeros(nfreqeval, n_oper, n_ener);
% for i = 1:n_ener
%    for j = 1:n_oper
%       for k = 1:nfreqeval
%          if(ifreq(i,(j-1)*nfreqeval+k)~= 0)
%             sres(k,j,i) = sres(k,j,i) - occ_sign(i,(j-1)*nfreqeval+k) * psiphi(:,(i-1)*n_oper+j)' * ...
%                (epsilon(:,:,ifreq(i,(j-1)*nfreqeval+k)) * fact1(i,(j-1)*nfreqeval+k) + ...
%                epsilon(:,:,ifreq(i,(j-1)*nfreqeval+k)+1) * fact2(i,(j-1)*nfreqeval+k))* ...
%                Dcoul * psiphi(:,(i-1)*n_oper+j) / vol;
%          end
%       end
%    end
% end
% disp('Real part done, go imag part.')
% %sres
% 
% %Integral contrib in imag frequency
% sw_imag_freq = zeros(nfreqimag, n_oper, n_ener);
% for i = 1:n_ener
%    for j = 1:n_oper
%       for k = 1:nfreqimag
%          sw_imag_freq(k,j,i) = psiphi(:,(i-1)*n_oper+j)' * epsilon(:,:,k+nfreqreal) * ...
%             Dcoul * psiphi(:,(i-1)*n_oper+j) / vol;
%       end
%    end
% end
% %sw_imag_freq
% sint = zeros(nfreqeval, n_oper, n_ener);
% imag_freq = imag(imag_freq) * ry2ev;
% switch options.cd_int_method
%    case 0
%       for i = 1:n_ener
%          for j = 1:n_oper
%             for k = 1:nfreqimag
%                if(k == 1)
%                   freqstart = imag_freq(k);
%                   freqend = (imag_freq(k) + imag_freq(k+1)) * 0.5;
%                   sint(:,j,i) = sint(:,j,i) + (sw_imag_freq(k,j,i) * ...
%                      (atan( freqend ./ wx(i,(j-1)*nfreqeval+1:j*nfreqeval)) - ...
%                      atan(freqstart ./ wx(i,(j-1)*nfreqeval+1:j*nfreqeval)))).';
%                elseif(k >=2 & k < nfreqimag)
%                   freqstart = (imag_freq(k-1) + imag_freq(k)) * 0.5;
%                   freqend = (imag_freq(k) + imag_freq(k+1)) * 0.5;
%                   sint(:,j,i) = sint(:,j,i) + (sw_imag_freq(k,j,i) * ...
%                      (atan( freqend ./ wx(i,(j-1)*nfreqeval+1:j*nfreqeval)) - ...
%                      atan(freqstart ./ wx(i,(j-1)*nfreqeval+1:j*nfreqeval)))).';
%                else
%                   freqstart = (imag_freq(k-1) + imag_freq(k)) * 0.5;
%                   freqend = (-imag_freq(k-1) + 3.0 * imag_freq(k)) * 0.5;
%                   sint(:,j,i) = sint(:,j,i) + (sw_imag_freq(k,j,i) * ...
%                      (atan( freqend ./ wx(i,(j-1)*nfreqeval+1:j*nfreqeval)) - ...
%                      atan(freqstart ./ wx(i,(j-1)*nfreqeval+1:j*nfreqeval)))).';
%                end
%             end
%          end
%       end
%    case 2
% 
%    case 3
% 
%    otherwise
%       error('Error: Invalid integration method')
% end
% 
% asxDtmp = zeros(nfreqeval,n_ener);
% achDtmp = zeros(nfreqeval,n_ener);
% achtDtmp = zeros(nfreqeval,n_oper,n_ener);
% %for i = 1:n_oper
%    for j = 1:n_ener
%       asxDtmp(:,j) = sum(sres(:,:,j), 2);
%       achDtmp(:,j) = sum(sint(:,:,j), 2)/pi;
%    end
% %end
% 
% %这部分不需要做，achtDtmp应为sint((nfreqeval+1)/2,1:n_oper,1:n_ener)
% if(0)
%    tmp = 1:1:nfreqeval;
%    tmp = abs((tmp-1) * options.delta_freq_step - options.max_freq_eval + TOL_SMALL);
%    for i = 1:n_ener
%       for j = 1:n_oper
%          diffmin = 1E12;
%          for k = 1:nfreqeval
%             if(tmp(k)<diffmin)
%                diffmin = tmp(k);
%                iw = k;
%             end
%          end
%          achtDtmp(iw,j,i) = sint(iw,j,i)/pi;
%       end
%    end
% else
%    achtDtmp((nfreqeval+1)/2,:,:) = sint((nfreqeval+1)/2,:,:)/pi;
% end
% %achtDtmp
% %asxDtmp
% %achDtmp
% 
% %asxtDyn = zeros(nfreqeval,n_ener);
% %achtDyn = zeros(nfreqeval,n_ener);
% %asxtDyn = sum(asxDtmp,2);
% %achtDyn = sum(achDtmp,2);
% 
% %asxtDyn
% %achtDyn
% Ex = zeros(n_ener,n_ener);
% 
% psioccnr = prod_states_gw(psir(:,1:nv),conj(psir(:,nv-nv_ener+1:nv+nc_ener)));
% psioccn = F*psioccnr;
% for iocc = nv-nv_oper+1:nv_oper
%   Ex = Ex + (psioccn(:,(iocc-1)*n_ener+1:iocc*n_ener))'*(Dcoul*psioccn(:,(iocc-1)*n_ener+1:iocc*n_ener));
% end
% 
% timeforSigma = toc(startSigma)
% 
% startEqp = tic;
% 
% fprintf('Unsymmetrized values')
% fac = 1.0/vol;
% E0=ev(nv-nv_ener+1:nv+nc_ener)*ry2ev
% Ex = real(diag(0.0 - fac * Ex)) * ry2ev
% Vxc = Vxc(nv-nv_ener+1:nv+nc_ener)
% Esx_x = (asxDtmp((nfreqeval+1)/2,:) * ry2ev)';
% Ech = (achDtmp((nfreqeval+1)/2,:) * ry2ev)';
% Sig = Ex + Esx_x + Ech
% Eqp0 = E0 - Vxc + Sig;
% save('Energy_ff.mat', 'E*', 'Sig');
% %E0=ev(nv-nv_ener+1:nv+nc_ener)
% %Ex = diag(real(0.0 - fac * Ex))
% %vxc = Vxc(nv-nv_ener+1:nv+nc_ener)/ry2ev
% %Esx_x = (real(asxDtmp((nfreqeval+1)/2,:)))'
% %Ech = (real(achDtmp((nfreqeval+1)/2,:)))'
% 
% %Computes and symmetrizes the quasiparticle spectrum
% dek = ev(nv-nv_ener+2:nv+nc_ener) - ev(nv-nv_ener+1:nv+nc_ener-1);
% tmp = find(dek<TOL_SMALL);
% if(size(tmp,1) == 0)
%    n = n_ener;
%    ndeg = 1:n_ener;
% elseif(size(tmp,1) == 1)
%    ndeg = zeros(1,n_ener-1);
%    n = n_ener-size(tmp,1);
%    if(tmp(1)~=1)
%       ndeg(1) = 2;
%       ndeg(2:-1) = 1;
%    else
%       ndeg(1:tmp(1)-1) = 1;
%       ndeg(tmp(1)) = 2;
%       ndeg(tmp(1)+1:-1) = 1;
%    end
% else
%    n = n_ener-size(tmp,1)
%    ndeg = zeros(1,n_ener-size(tmp,1));
%    error('Error: Not support Now!')
% end      
% 
% asigt_cor = achDtmp*ry2ev;
% asigt2 = asxDtmp*ry2ev;
% asigt_corb = asigt_cor + asigt2;
% asigt = asigt_corb + Ex.';
% asig = (asigt2((nfreqeval+1)/2,:) + asigt_cor((nfreqeval+1)/2,:)).' + Ex;
% %achDtmp*ry2ev
% 
% efsto = E0 - Vxc + asig;
% %freq0
% freqs = zeros(nfreqeval, n_ener);
% iw = (1:nfreqeval).';
% for i = 1:n_ener
%    freqs(:,i) = freq0(i) + (iw-1)*options.delta_freq_step;
% end
% 
% fmin = zeros(nfreqeval, n_ener);
% for i = 1:n_ener
%    fmin(:,i) = E0(i) + asigt(:,i) - Vxc(i) -freqs(:,i);
% end
% 
% Eqp1 = zeros(n_ener,1);
% solns = zeros(nfreqeval,n_ener);
% neqp1 = zeros(n_ener,1);
% if(nfreqeval < 2)
%    Eqp1 = eqp0;
%    neqp1 = 0;
% else
%    nsols = 0;
%    isol = 1;
%    dw = freqs(2,:) - freqs(1,:);
%    for i = 1:n_ener
%       for j = 2:nfreqeval
%          if(real(fmin(j-1,i))>0 & real(fmin(j,i))<0)
%             nsols = nsols + 1;
%             solns(nsols, i) = freqs(j-1,i) + (fmin(j-1,i) * dw(i)) / (real(fmin(j-1,i)) - real(fmin(j,i)));
%             if(abs(real(solns(nsols)-eqp0)) < abs(real(solns(isol)-eqp0)) | nsols ==1)
%                isol = nsols;
%             end
%          end
%       end
%       neqp1(i) = nsols;
%       if(nsols == 0)
%          if(real(fmin(1,i))<0 & real(fmin(1,i))>real(fmin(2,i)))
%             solns(1,i) = (fmin(2,i)*freqs(1,i) - fmin(1,i)*freqs(2,i)) / (real(fmin(2,i)) - real(fmin(1,i)));
%             neqp1(i) = -1;
%          elseif(real(fmin(nfreqeval,i))>0 & real(fmin(nfreqeval-1,i))>real(fmin(nfreqeval,i)))
%             solns(1,i) = (fmin(nfreqeval-1,i)*freqs(nfreqeval,i) - fmin(nfreqeval,i)*freqs(nfreqeval-1,i)) / (real(fmin(nfreqeval-1,i)) - real(fmin(nfreqeval,i)));
%             neqp1(i) = -2;
%          else
%             solns(nsols+1,i) =  Eqp0(i);
%             neqp1(i) = 0;
%          end
%       end
%       Eqp1(i) = solns(isol,i);
%    end
% end
% 
% fprintf('Symmetrized values from band-averaging');
% 
% Cor = (asigt_corb((nfreqeval+1)/2,:)).'
% Eqp0
% Eqp1
% 
% timeforEqp = toc(startEqp)
