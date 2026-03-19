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

% 
tmp = GWinfor.tmp_devel;
qindx_S = tmp.qindx_S;
nbz = tmp.nbz;
nibz = tmp.nibz;
kpt = tmp.kpt_Cart;
kptbz = tmp.kptbz_Cart;
sstar = [tmp.bz2ibz, tmp.bz2rot];




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
for ik = 1:nibz
  for ib = nbmin:nbmax
    for ispin = 1:nspin
      for iq = 1:nbz
        iqibz = sstar(iq, 1);
        iqs = sstar(iq, 2);
                
        ikp_bz = qindx_S(ik, iq, 1);
        is = qindx_S(ik, iq, 2);
        ikp_ibz = sstar(ikp_bz, 1);
        ikp_rot = sstar(ikp_bz, 2);
        
        vcoul_q = coulomb_main(GWinfor, config, iqibz);
        % vcoul_q(abs(vcoul_q) < 1e-10) = GWinfor.coulG0;
        % error("under construction")
        for ob = 1:size(GWinfor.occupation, 1)
          occ = GWinfor.occupation(ob, ikp_ibz, ispin);
          if (occ < 1e-6)
            continue;
          end

          param = [];
          param.is = [ib, ik, 1, ispin];
          param.os = [ob, ikp_ibz, ikp_rot, ispin];
          param.qs = [is, iqibz, iqs];

          if (ob == 1)
              disp(param)
              %
              kleft = tmp.kpt_RLU(ik, :);
              kright_bz = tmp.kptbz_RLU(ikp_bz, :);
              kright_ibz = tmp.kpt_RLU(ikp_ibz, :);
              Sk = tmp.rot_mtrx_RLU_G{ikp_rot};
              krightSk = kright_bz * Sk;
              qbz = tmp.kptbz_RLU(iq, :);
              qibz = tmp.kpt_RLU(iqibz, :);
              Sq = tmp.rot_mtrx_RLU_G{iqs};
              iGo = tmp.Ggrid_RLU(is, :);
              out1 = kleft - kright_bz - qibz*Sq; % k - kp - S*q_ibz
              out2 = iGo;
              % disp( out1 )
              % disp( out2 ) 
;
          end

          ngrho_left = mtxel(GWinfor, param);

          % call mtxel to get <..|..|..>
          DL_vol = det(GWinfor.supercell);
          RL_vol = (2*pi)^3 / DL_vol;
          d3q_factor = RL_vol / nbz;
          q_weight = d3q_factor / (2*pi)^3;
          Ry2eV = 13.6059;
          % tmp = sum(occ .* vcoul_q .* abs(ngrho_left).^2) / vol  * q_weight / 2;
          coeff = 4*pi; %  --> Ha
          % vcoul_q_tmp = vcoul_q / (8*pi) * q_weight;
          Ex_t = sum(vcoul_q .* abs(ngrho_left).^2);
          fprintf('nb = %3d, out = %.3e\n', ob, Ex_t / 8 / pi)
          Ex(ib, ik, ispin) = Ex(ib, ik, ispin) + Ex_t;
          % call matrix multiplication
          % vcoul_q_tmp(1:20)
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