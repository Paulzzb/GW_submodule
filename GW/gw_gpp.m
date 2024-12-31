function gw_gpp(ksinfo, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is used to calculate GW quasiparticle energies under GPP approximatio.
%
% Input:
% ksinfo: A structure that contains ground state calculation results
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

testInput = true;
testInput = false;
testflag1 = true;
testflag1 = false;
testflag2 = true;
testflag2 = false;
testflag3 = true;

flagscalarize = false;
isCPLX = true;

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
    strEval = sprintf('%s = %.16f;', fieldname, value)
    eval(strEval);
	end
%    fprintf('Field: %s, Value: %s\n', fieldname, num2str(value));
end

if testInput
  ksinfor2 = load("Siksinfo.mat");
  ksinfor2 = ksinfor2.ksinfo;
  ksinfo.Z = ksinfor2.Z;
  ksinfo.ev = ksinfor2.ev;
  % ksinfo.aqs = ksinfor2.aqs;
end

% dcoul = ksinfo.coulG(:, 4);

% limitone =  1.0 / (4.0 * TOL_SMALL);
% limittwo = 0.50.^2;
% if isfield(options, 'gpp_brodening')
%   limittwo = gpp_brodening .^ 2;
% end
% sexcut = 4.0;
  
ksinfo.Z = ksinfo.Z * sqrt(ksinfo.vol); % Turn to \int_V \abs{psi(r)}^2 dr = 1.
Z     = ksinfo.Z;
F     = ksinfo.F;
ev    = ksinfo.ev;
Vxc   = ksinfo.Vxc;
dcoul = ksinfo.coulG(:,4);
dcoul(1) = ksinfo.coulG0;
Dcoul = spdiags(ksinfo.coulG(:,4), 0, ng, ng);
Dcoul(1, 1) = ksinfo.coulG0;
gvec = ksinfo.gvec;
if isempty(ksinfo.aqs)
  aqsFlag = false;
else
  aqsFlag = true;
end
% Dcoul(1,1) = ksinfo.coulG0;
% fprintf('nr = %d, ng = %d, n_oper = %d\n', nr, ng, n_oper);



% Normalize the wavefunc in Fourier space.
% startZ = tic;
% timeforZ = toc(startZ);
% startC2R = tic;
% psir = zeros(nr, nv+nc);
% for iband = 1:nv+nc
%   fftbox1 = put_into_fftbox(Z(:, iband), gvec.idxnz, gvec.fftgrid);
%   fftbox1 = gvec.nfftgridpts / ksinfo.vol * do_FFT(fftbox1, gvec.fftgrid, 1);
%   psir(:, iband) = reshape(fftbox1, gvec.nfftgridpts, []);
% end
% timeforC2R = toc(startC2R);
% fprintf('Time for psir = %.4f.\n', timeforC2R)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start calculation

startforW = tic;



Mgvc = zeros(ng, nc_oper);
Eden = (1) ./ (kron(ev(nv-nv_oper+1:nv), ones(nc_oper,1)) ... 
            - kron(ones(nv_oper, 1), ev(nv+1:nv+nc_oper)));
inveps = zeros(ng, ng);
scal = 4.0;
% scal = 4.0 / ksinfo.vol; % for kpoints & spin, change it according to
            % epsilon_main.f90 and chi_summation.f90
for ind_nv = nv-nv_oper+1:nv
  if aqsFlag
    Mgvc = ksinfo.aqs{ind_nv}(:, nv+1:nv+nc_oper);
  else    
    Mgvc = mtxel_sigma(ind_nv, ksinfo, options.Groundstate, (nv+1:nv+nc_oper));
  end
  Mgvc = conj(Mgvc);
  eden = 1 ./ (ev(ind_nv) - ev(nv+1:nv+nc_oper));
  inveps = inveps + scal * Mgvc * diag(eden) * Mgvc';
end
if testflag1
  fprintf("||chi||_F = %.3e.\n", norm(inveps, 'fro'))
  tocompare = readmatrix('chi.csv');
end
% inveps = eye(ng) - Dcoul * inveps;
inveps = eye(ng) - Dcoul * inveps / vol;
if testflag1
  fprintf("||epsilon||_F = %.3e.\n", norm(inveps, 'fro'))
  tocompare = readmatrix('epsilon.csv');
  fprintf("||epsilon||_F = %.3e.\n", norm(tocompare, 'fro'))
end
if testflag1
  dlmwrite('epsilon_.csv', inveps, 'precision', 16);
end
inveps = inv(inveps);
if testflag1
  fprintf("||epsilon^-1||_F = %.3e.\n", norm(inveps, 'fro'))
  tocompare = readmatrix('inveps.csv');
  fprintf("||epsilon^-1||_F = %.3e.\n", norm(tocompare, 'fro'))
end
if testflag1
  dlmwrite('inveps_.csv', inveps, 'precision', 16);
end
I_inveps = eye(ng) - inveps;
if testflag1
  fprintf("||I - epsilon^-1||_F = %.3e.\n", norm(I_inveps, 'fro'))
  tocompare = readmatrix('W.csv');
  fprintf("||I - epsilon^-1||_F = %.3e.\n", norm(tocompare, 'fro'))
end
if testflag1
  tocompare = readmatrix('W.csv');
  norm(tocompare - I_inveps)
  dlmwrite('W_.csv', inveps, 'precision', 16);
end

% Use result from old one to replace I_inveps
if testflag1
  I_inveps = tocompare;
end

TimeofW = toc(startforW);
fprintf("Time to get I - inv(epsilon) = %.3f.\n", TimeofW);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate wtilde
startwtilde = tic;
% flagscalarize = true;


for igcol = 1:ng
  Omega2 =  wpeff(ksinfo, options, igcol); % There should not be 1/vol here...
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
timewtilde = toc(startwtilde);
fprintf("Time for wtilde = %.3f.\n", timewtilde);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare wtilde_array
if testflag1
  wtilde_array_bgw = readmatrix('wtilde_array.csv');
  for igcol = 1:ng
    nonzero_wtilde_array = find(wtilde_array(:, igcol)); 
    nonzero_wtilde_array_bgw = find(wtilde_array_bgw(:, igcol)); 
    if length(nonzero_wtilde_array_bgw) ~= length(nonzero_wtilde_array);
      fprintf('Different number of non-zero indices!\n');
      fprintf('%d : bgw, %d : Ours\n', length(nonzero_wtilde_array_bgw), length(nonzero_wtilde_array));
      warning();
      continue
    end
    if any(nonzero_wtilde_array ~= nonzero_wtilde_array_bgw)
      disp(nonzero_wtilde_array)
      disp(nonzero_wtilde_array_bgw)
      error('Even the non-zero indices are not the same in wtilde_array');
    else
      wtilde_array_diff = 1./ wtilde_array(nonzero_wtilde_array, igcol) ...
                        - 1./wtilde_array_bgw(nonzero_wtilde_array_bgw, igcol);
      fprintf('igcol = %d, difference in wtilde_array = %f\n', ...
              igcol, norm(wtilde_array_diff));
      if (norm(wtilde_array_diff) / norm(1./wtilde_array(nonzero_wtilde_array, igcol))) >= 1e-5
        fprintf('Norm difference of 1/wtilde_array = %f.\n', ...
                (norm(wtilde_array_diff) / norm(1./wtilde_array(nonzero_wtilde_array, igcol))));
        warning('wtilde_array seems wrong!')
        pause(0.1)
      end;
    end  
  end  

  % tocompare = readmatrix('wtilde_array.csv');

  % out1 = norm(wtilde_array, 'fro');
  % out2 = norm(tocompare, 'fro');
  % fprintf("||wtilde||_F = %.3e, ||wtilde_real||_F = %.3e.\n", out1, out2);
  % if norm(wtilde_array - tocompare) > 1e-6
  %   fprintf("wtilde_array may wrong!\n");
  % end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Construct Sigma
startSelfEnergy = tic;
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
  if aqsFlag
    aqsntemp = ksinfo.aqs{ibandouter_aqsntemp};
  else
    aqsntemp = mtxel_sigma(ibandouter, ksinfo, options.Groundstate, ...
              (nv-nv_oper+1:nv+nc_oper));
  end
  % THIS LINE ONLY to CHECK correctness!!!
  % COEFFICIENTS should be CHNAGES later!!!
  ssxflag = false;
  schflag = false;
  schttflag = false;
  asxtemp = ZERO;
  achtemp = ZERO; 
  for ibandinner = nv-nv_oper+1 : nv+nc_oper % line 1214 of mtxel_cor
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
        
        % if flagscalarize
        %   wtilde = wtilde_array(:, igcol);
        %   wtilde2 = wtilde.^2;
        %   wdiff = ones(ng, 1) * wxt;
        %   cden = wdiff; rden = cden .* conj(cden); rden = 1.0 ./ rden;
        %   delw = wtilde .* conj(cden) .* rden;
        %   delwr = delw .* conj(delw);
        % end
        for igrow = 1:ng % line 1508
          wtilde = wtilde_array(igrow, igcol);
          wtilde2 = wtilde .^ 2;
          Omega2 = wtilde2 * I_inveps(igrow, igcol);
          wdiff = wxt - wtilde;
          cden = wdiff;
          rden = cden * conj(cden);
          rden = 1.0 / rden;
          delw = wtilde * conj(cden) * rden; % just help it to inverse a conmplex number.
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
        end % for igrow = 1:ng
        sch_array = sch_array + 0.5 * scht * conj(aqsntemp(igcol, ibandinner_aqsntemp));
      end % if flagocc
  
      if (flagocc == true)
        asxtemp = asxtemp - ssx_array * occ * dcoul(igcol) / ksinfo.vol;
      end
      achtemp = achtemp + sch_array * dcoul(igcol) / ksinfo.vol;
    end % for igcol = 1:ng
  
  end % for ibandinner
  
  asxt = asxt + asxtemp;% line 1280, mtxel_cor
  acht = acht + achtemp;
  Esx_x(ibandouter_aqsntemp) = asxt;
  Ech(ibandouter_aqsntemp)   = acht;
  Ex(ibandouter_aqsntemp) = sum(diag(aqsntemp(:, 1:nv_oper)' * Dcoul ...
                            * aqsntemp(:, 1:nv_oper) / vol));
  
end % for ibandouter

% for ioper = nv-nv_oper+1:nv
%   if aqsFlag
%     Mgvn = ksinfo.aqs{ioper}(:, nv-nv_ener+1:nv+nc_ener);
%   else
%     Mgvn = mtxel_sigma(ioper, ksinfo, options.Groundstate, (nv-nv_ener+1:nv+nc_ener));
%   end
%   Mgvn = conj(Mgvn);
%   Ex = Ex + Mgvn' * Dcoul * Mgvn / ksinfo.vol;
% end

TimeofSigma = toc(startSelfEnergy);


fprintf("Time to self-energies = %.3f.\n", TimeofSigma);


% Calculate energies
GWenergy.ev = ev(nv-nv_ener+1:nv+nc_ener) * ry2ev;
GWenergy.Ex = real((0.0 + Ex) * ry2ev);
GWenergy.Esx_x = real((0.0 + Esx_x) * ry2ev);
GWenergy.Ech = real((1.0 * Ech * ry2ev));
GWenergy.Sig = real(GWenergy.Ex + GWenergy.Esx_x + GWenergy.Ech);
GWenergy.Vxc = Vxc(nv-nv_ener+1:nv+nc_ener) * ry2ev;
GWenergy.eqp = GWenergy.ev - GWenergy.Vxc + GWenergy.Sig;




timeforGW = toc(startGW);
save(options.GWCal.fileName, 'GWenergy', 'time*') 
fprintf('Time for GW Calculation under COHSEX approximation = %f.\n', timeforGW);
fprintf('Result is saved in %s.\n', options.GWCal.fileName)
