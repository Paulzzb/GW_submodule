function Ex = gw_x_k(GWinfor, config)


msg = sprintf('[Exchange] Start computing Σ_x (exchange part)...\n');
QPlog(msg, 0);
tStart = tic;

default_Constant = constant_map();
nameConstants = fieldnames(default_Constant);
for i = 1:numel(nameConstants)
  eval(sprintf('%s = %.16f;', nameConstants{i}, default_Constant.(nameConstants{i})));
end



% Basic setup
ng = GWinfor.gvec.ng; % corresponds to Dcoul
nbmin = config.SYSTEM.energy_band_index_min;
nbmax = config.SYSTEM.energy_band_index_max;
nv = find(GWinfor.occupation > 1 - TOL_SMALL, 1, 'last');
nspin = GWinfor.nspin;
nspinor = GWinfor.nspinor;

vol = GWinfor.vol;
gvecCoul = GWinfor.gvec;
psir = GWinfor.psir;

symminfo = GWinfor.symminfo;
bz_samp = GWinfor.bz_samp;

nbz = bz_samp.nbz;
kpt = bz_samp.kpt;
nibz = bz_samp.nibz;
kptbz = bz_samp.kptbz;






msg = sprintf('[Exchange] Start computing Σ_x (exchange part)...\n');
QPlog(msg, 0);
tStart = tic;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Main part Ex with k-points

msg = sprintf('[Exchange] Using standard Σ_x calculation.\n');
QPlog(msg, 0);

% error('GW_x_k.m is not ready yet.');

Ex = zeros(nbmax-nbmin+1, nibz, nspin);
tStandard = tic;
for ib = nbmin:nbmax
  for ik = 1:nibz
    for ispin = 1:nspin
      for iq = 1:nbz
        iq_ibz = bz_samp.kbz2kibz_ind_kbz(iq);
        iqs = bz_samp.kbz2kibz_ind_rotation(iq);
        vcoul_q = GWinfor.coulG_list;

        ikp_bz = bz_samp.qindx_S(ik, iq, 1);
        is = bz_samp.qindx_S(ik, iq, 2);
        ikp_ibz = bz_samp.kbz2kibz_ind_kbz(ikp_bz);
        ikp_rot = bz_samp.kbz2kibz_ind_rotation(ikp_bz);
        iGo = bz_samp.iGolist(ikp_ibz);
        
        % error("under construction")
        for ob = 1:size(GWinfor.occupation, 1)
          param = [];
          param.is = [ib, ik, 1, ispin];
          param.os = [ob, ikp_ibz, ikp_rot, ispin];
          param.qs = [is, iq_ibz, iqs];

          ngrho_left = mtxel(GWinfor, param);
          occ = GWinfor.occupation(ob, ikp_ibz, ispin);

          % call mtxel to get <..|..|..>
          Ex(ib, ik, ispin) = Ex(ib, ik, ispin) + sum(occ .* vcoul_q{iq_ibz} .* abs(ngrho_left).^2) / vol / nbz ;
          % call matrix multiplication
        end
      end
    end % ispin
  end % ik
end % ib
msg = sprintf('[Exchange] Standard loop completed in %.2f seconds.\n', toc(tStandard));
QPlog(msg, 1);

%% End main part
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ex = - real(diag(Ex));
msg = sprintf('[Exchange] Finished. Total time: %.2f seconds.\n', toc(tStart));
QPlog(msg, 0);

end % EOF