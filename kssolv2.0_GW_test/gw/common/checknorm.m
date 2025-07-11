function checknorm(wfn, ik, ispin)
norm = vecnorm(wfn);   %波函数的结构是options.X0.wavefuncell{k, spin}.psi(NG,iband)
norm_residual = norm - 1;
if all(norm_residual(:) > 1000*eps)
    iband = find(norm_residual(:) > 1000*eps)
    error('norm of wavefunction iband ik=%d ispin=%d is not 1',ik, ispin)
end
end