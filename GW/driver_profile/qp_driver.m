function qp_driver(input_dir)
% qp_driver -> driver to perform quasi-particle calculation

def = filename_map();
fName = fullfile(input_dir, def.GWinput);
TEMP = load(fName);

GWinfo = TEMP.GWgroundstate;
config = TEMP.config;

init_log(config.CONTROL.log_level);  % initialize log
GWlog('[QPdriver] QP driver started', 0);
startQP = tic;

% Prepare real-space wavefunction from psig
GWlog('[QPdriver] Converting wavefunction from reciprocial space to real space ...', 1);
GWinfo.psir = get_wavefunc_real(GWinfo.psig, GWinfo.Ggrid4psig);
GWlog('[QPdriver] Wavefunction in real space prepared.', 2);

% Initialize energy structure
GWenergy = QPenergy(GWinfo, config);
GWlog('[QPdriver] Quasiparticle energy structure initialized.', 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GW method start
% Begin main GW calculation
if config.CONTROL.isgw
  GWenergy = qpgw(GWinfo, config);
end % config.&CONTROL.isgw

% GW method done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 3: Energy shift to account for degeneracy, etc.
GWlog('[QP driver] Post-processing QP energy shift...', 1);
GWenergy = shiftenergy(GWenergy);
GWlog('[QP driver] Energy shift completed.', 2);

% Step 4: Compute final E_qp and output
GWlog('[QP driver] Computing final quasiparticle energies...', 1);
GWenergy = getEqp(GWenergy);

GWlog('[QP driver] Saving QP results to output...', 1);
GWfout(GWenergy);

% Finish timing and wrap up
timeQP = toc(startQP);
msg = sprintf('[QP driver] Quasiparticle calculation finished. Total time: %.2f seconds.', timeQP);
GWlog(msg, 0);


% Output
GWfout(GWenergy);



end % main function

