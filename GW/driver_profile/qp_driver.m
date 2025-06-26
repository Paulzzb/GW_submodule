function qp_driver(input_dir)
% qp_driver -> driver to perform quasi-particle calculation
def = filename_map();
fName = fullfile(input_dir, def.GWinput);
TEMP = load(fName);

GWinfo = TEMP.GWgroundstate;
config = TEMP.config;


% Initialize Log information
cleanup = QPlog_push('qp_driver');
QPlog_showtag(true);
QPlog_verbose(config.CONTROL.log_level);
QPlog('QP driver started', 0);
startQP = tic;




% Prepare real-space wavefunction from psig
QPlog('Converting wavefunction from reciprocial space to real space ...', 1);
GWinfo.psir = get_wavefunc_real(GWinfo.psig, GWinfo.Ggrid4psig);
QPlog('Wavefunction in real space prepared.', 2);

% Initialize energy structure
GWenergy = QPenergy(GWinfo, config);
QPlog('Quasiparticle energy structure initialized.', 2);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GW method start
% Begin main GW calculation
if config.CONTROL.isgw
  GWenergy = qpgw(GWinfo, config);
end % config.&CONTROL.isgw

% GW method done
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step 3: Energy shift to account for degeneracy, etc.
QPlog('Post-processing QP energy shift...', 1);
GWenergy = shiftenergy(GWenergy);
QPlog('Energy shift completed.', 2);

% Step 4: Compute final E_qp and output
QPlog('Computing final quasiparticle energies...', 1);
GWenergy = getEqp(GWenergy);

msg = sprintf('Saving QP-energies results to output file %s...', ...
              GWenergy.fout);
QPlog(msg, 1);
GWfout(GWenergy);

% Finish timing and wrap up
timeQP = toc(startQP);
msg = sprintf('Quasiparticle calculation finished. Total time: %.2f seconds.', timeQP);
QPlog(msg, 0);



end % main function

