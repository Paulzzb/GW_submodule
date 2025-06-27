function qp_postprocess(GWenergy)
  cleanup = QPlog_push('qp-postprocess');

  QPlog('begin post-processing of QP energy...', 0);

  % step 3: energy shift to account for degeneracy, etc.
  QPlog('post-processing QP energy shift...', 1);
  GWenergy = shiftenergy(GWenergy);
  QPlog('energy shift completed.', 2);
  
  % step 4: compute final e_QP and output
  QPlog('computing final quasiparticle energies...', 1);
  GWenergy = getEqp(GWenergy);
  
  msg = sprintf('saving QP-energies results to output file %s...', ...
                 GWenergy.fout);
  QPlog(msg, 1);
  GWfout(GWenergy);
  
  QPlog('QP post-processing completed.', 0);



end
