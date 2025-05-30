function nm_Xomega_nm = fourcenterintegral(GWinfo, options, Wflag, ...
                      n_start_end, m_start_end, omega_list, varargin)

% Calculate nm_X_nm(nid, mid) = <nm | W(q=0;omega)| nm >
% for all n in n_start_end, m in m_start_end
% If pattern is provided

if (nargin < 6)
  error("First 6-th arguments are required");
end

% Some controls

isisdf = options.ISDFCauchy.isISDF;

eta = 0.0;

% Allocate space
nstart = n_start_end(1);
nend   = n_start_end(2);
mstart = m_start_end(1);
mend   = m_start_end(2);
Nn = nend - nstart + 1;
Nm = mend - mstart + 1;
mlist = mstart:mend;
Nw = length(omega_list);
if Wflag == 0
  Nw = 1;
end
nm_Xomega_nm = zeros(Nn, Nm, Nw);


patternflag = 0;

for i = 1:2:length(varargin)
    switch lower(varargin{i})
        case 'pattern'
            pattern = varargin{i+1};
            patternflag = 1; % Provided pattern
    end
end

if (patternflag == 1)
  if ~isequal(size(pattern), [Nn, Nm, Nw])
    fprintf('Pattern size is not consistent with n_start_end, m_start_end, and omega_list!')
    fprintf("In case Wflag = 0, size of pattern is [Nn, Nm, 1]");
    error("Error in fourcenterintegral.m");
  end
else
  pattern = ones(Nn, Nm, Nw); % Default pattern
end




% 0. Initialize
nameConstants = fieldnames(options.Constant);
for i = 1:numel(nameConstants)
  fieldname = nameConstants{i};
  value = options.Constant.(fieldname);    
  if ~isempty(value)
    strEval = sprintf('%s = %.16f;', fieldname, value);
    eval(strEval);
  end
end


% Initialize other values
% before it is 1-norm based, now we need \sqrt(vol) based.
% GWinfo.Z = GWinfo.Z * ry2ev;
mi = sqrt(-1);
ev    = GWinfo.ev;
gvec  = GWinfo.gvec;
Dcoul = spdiags(GWinfo.coulG(:,4), 0, ng, ng);
Dcoul(1, 1) = GWinfo.coulG0;

% Energies use unit 'ev' in this code.
GWinfo.Z = GWinfo.Z * sqrt(vol);
ev = ev * ry2ev;
Dcoul = Dcoul * ry2ev;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% No ISDF
if ~isisdf
  for ifreq = 1:Nw
    omega = omega_list(ifreq);
    if (abs(real(omega)) < 1e-8)
      ishermW = 1;
    else
      ishermW = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if Wflag
      % No ISDF, use the screened Coulomb matrix
      % Calculate W
      epsilon = zeros(ng, ng);
      for ind_nv = nv-nv_oper+1:nv
        Mgvc = mtxel_sigma(ind_nv, GWinfo, options.Groundstate, ...
                (nv+1:nv+nc_oper));
        Mgvc = conj(Mgvc);
        
        Eden = ev(ind_nv) - ev(nv+1:nv+nc_oper);
        edenDRtmp = (-1.0 ./ (omega - Eden - mi*eta) ...
        + 1.0 ./ (omega + Eden + mi*eta));
        epsilon = epsilon + 2*Mgvc*(edenDRtmp.*Mgvc') / vol;
      end % for ind_nv
      if ishermW
        epsilon = (Dcoul/vol)^(-1) - epsilon*vol;
        % Here epsilon is supposed to be an Hermite matrix.
        epsilon = tril(epsilon, -1) + tril(epsilon, -1)' + diag(real(diag(epsilon)));
        [L, D] = ldl(epsilon); rsqrtD = diag(sqrt(diag(D)).^(-1));
      else
        epsilon = eye(ng) - Dcoul * epsilon;
        epsilon = inv(epsilon);
        W = (eye(ng) - epsilon) * Dcoul / vol;
        Wnorm = norm(W, 'fro') ;
        fprintf(' Wnorm = %12.6f\n',  Wnorm);
      end
      
      for n = nstart:nend
        % Use pattern to check, for a given n, mlist <nm| 
        tmp = pattern(n-nstart+1, :, ifreq);
        indm = find(tmp);
        mlisttmp = mlist(indm);
        Mgvc = mtxel_sigma(n, GWinfo, options.Groundstate, ...
                          mlisttmp);
        Mgvc = conj(Mgvc);
        if ishermW
          out_list = sum(Dcoul/vol*abs(Mgvc).^2)';
          Mgvc = rsqrtD * (L\Mgvc); 
          out_list = out_list - sum(abs(Mgvc).^2)';
          nm_Xomega_nm(n-nstart+1, indm, ifreq) = out_list;
        else
          for i = 1:length(mlisttmp)
            nm_Xomega_nm(n-nstart+1, indm(i), ifreq) = Mgvc(:, i)'*W*Mgvc(:, i);
          end
        end
      end
    else % Wflag
    % No ISDF, just use the bare coulomb matrix
      for n = nstart:nend
        % Use pattern to check, for a given n, mlist <nm| 
        tmp = pattern(n-nstart+1, :, 1);
        indm = find(tmp);
        mlisttmp = mlist(indm);
        Mgvc = mtxel_sigma(ibandener, GWinfo, options.Groundstate, ...
                          mlisttmp);
        Mgvc = conj(Mgvc);
        out_list = sum(Dcoul/vol*abs(Mgvc).^2)';
        nm_Xomega_nm(n-nstart+1, indm) = out_list;
      end
    end % Wflag 
  end % for ifreq
end % ~isisdf
% Finished

if isisdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
  
  if Wflag
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
    timeforISDF = toc(startISDF);
  
    % Finish ISDF
  end

  % Prepare <helper|V|helper>
  if Wflag
    vcVvc = vcgzeta_mu' * Dcoul * vcgzeta_mu / vol;
    vcVnn = vcgzeta_mu' * Dcoul * ssgzeta_mu / vol;
  else
    vnVvn = vsgzeta_mu' * Dcoul * vsgzeta_mu / vol;
    U = chol(vnVvn);
  end
end % End of ISDF


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use ISDF
if isisdf
  if Wflag
    for ifreq = 1:Nw
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ISDF, use the screened Coulomb matrix
    % Prepare W(q; omega)
      omega = omega_list(ifreq); 
      if (abs(real(omega)) < 1e-8)
        ishermW = 1;
      else
        ishermW = 0;
      end
      epsKernel = zeros(vcrank_mu, vcrank_mu);
      for ind_nv = nv-nv_oper+1:nv
        Mgvc = (psir(vcind_mu, ind_nv)) .* conj(psir(vcind_mu, nv+1:nv+nc_oper)); 
        Eden = ev(ind_nv) - ev(nv+1:nv+nc_oper);
        edenDRtmp = (-1.0 ./ (omega - Eden - mi*eta) ...
        + 1.0 ./ (omega + Eden + mi*eta));
        epsKernel = epsKernel + Mgvc*diag(edenDRtmp)*Mgvc';
      end % for ind_nv     
      epsKernel = inv(epsKernel)/2 - vcVvc; 

      if ishermW
        % Here epsKernel is supposed to be an Hermite matrix.
        epsKernel = tril(epsKernel, -1) + tril(epsKernel, -1)' + diag(real(diag(epsKernel)));
        [LKernel, DKernel] = ldl(epsKernel);
        dKernel = diag(DKernel); 
        right_nnWnn = LKernel \ vcVnn;
      else
        WKernel = -inv(epsKernel);
        clear epsKernel
      end %ishermW

      for n = nstart:nend
        % Use pattern to check, for a given n, mlist <nm| 
        tmp = pattern(n-nstart+1, :, ifreq);
        indm = find(tmp);
        if isempty(indm)
          continue;
        end
        mlisttmp = mlist(indm);
        Mgvc = psir(ssind_mu, n) .* conj(psir(ssind_mu, mlisttmp));
        if ishermW
          Mgvc = right_nnWnn*Mgvc; 
          out_list = - sum(dKernel.^(-1) .* abs(Mgvc).^2);
          nm_Xomega_nm(n-nstart+1, indm, ifreq) = out_list;
        else
          Mgvc = vcVnn * Mgvc;
          for i = 1:length(mlisttmp)
            nm_Xomega_nm(n-nstart+1, indm(i), ifreq) = Mgvc(:, i)'*WKernel*Mgvc(:, i);
          end
        end
      end
    end %for ifreq
  else
    % ISDF, bare Coulomb matrix 
    for n = nstart:nend
      tmp = pattern(n-nstart+1, :, 1);
      indm = find(tmp);
      mlisttmp = mlist(indm);
      Mgvn = (psir(vsind_mu, n)) .* conj(psir(vsind_mu, mlisttmp));
      Mgvn = conj(Mgvn);
      out_list = sum(U \ abs(Mgvn).^2)';
      nm_Xomega_nm(n-nstart+1, indm) = out_list;
    end
  end % if Wflag
end % isisdf

end % function