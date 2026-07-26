function qp_driver_package(stageFilePath)
% qp_driver -> driver to perform quasi-particle calculation

service_reset_persistent();
packages_reset_persistent();

def = filename_map();
test_stage(stageFilePath);


% Initialize Log information
cleanup = QPlog_push('qp_driver');
QPlog_showtag(true);
QPlog_verbose(1);
QPlog('QP driver package started', 0);
startQP = tic;




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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [Developer Hook] Insert your custom module calls below
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example: call your module if enabled in config
if isfield(config.CONTROL, 'enable_your_module') && config.CONTROL.enable_your_module
  QPlog('Your module is enabled. Starting execution...', 1);
  GWenergy = your_kernel(GWinfo, config);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% finish timing and wrap up
timeQP = toc(startQP);
msg = sprintf('quasiparticle calculation finished. total time: %.2f seconds.', timeQP);
QPlog(msg, 0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Post-processing
QP_postprocess(GWenergy)
% % step 3: energy shift to account for degeneracy, etc.
% QPlog('post-processing QP energy shift...', 1);
% GWenergy = shiftenergy(GWenergy);
% QPlog('energy shift completed.', 2);

% % step 4: compute final e_QP and output
% QPlog('computing final quasiparticle energies...', 1);
% GWenergy = getEqp(GWenergy);

% msg = sprintf('saving QP-energies results to output file %s...', ...
%               GWenergy.fout);
% QPlog(msg, 1);
% GWfout(GWenergy);

% % finish timing and wrap up
% timeQP = toc(startQP);
% msg = sprintf('quasiparticle calculation finished. total time: %.2f seconds.', timeQP);
% QPlog(msg, 0);


end % main function

