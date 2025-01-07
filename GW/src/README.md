# Structure in source code
```


├── BSE: Not used currently !!!

├── Common: Contains commonly used functions and classes.
|   The following are classes in GW module.
|   ├── options: options for GW module, @GWOptions
|   The following are widely used functions.
|   ├── do_FFT & get_from_fftbox & put_into_fftbox:
|       a set of functions to do FFT.
|   ├── gvec_to_fft_index: N x N x N --> N^3, in kssolv manner.
|                          NOT USED NOW.
|   ├── mtxel_sigma: Calculate <m|e^{ig\hat{r}}|n>.
|   ├── find_vector: find index iout to satisfy gvec(iout, :) = kk.
|                    Used in exact_CH and gw_gpp.
|   The following are used in fullfrequency calculation.
|   ├── freqgen: Generate frequency grid and corresponding weight function.
|                Used in gw_fullfreq_cd*.
|   └── GaussLegendre: Gauss-Legendre quadrature implementation.
|                Used in gw_fullfreq_cd_int*.


├── ISDF: Suppose that we are calculating
    Phi_i(r)*Psi_j(r) \approx \sum_{mu} p_mu(r) phi_i(r_mu)*psi_j(r_mu)
    p_mu(g) = \int_{R^3} dr exp^(i g r) p_mu(r).
|   ├── isdf_main: Main function to generate ISDF
|   ├── isdf_indices: Generate corresponding indices for r_mu, phi_i(r_mu), psi_j(r_mu).
|   ├── k_means: An implementation of k-means clustering.
|   ├── k_means_hpc: An implementation of k-means clustering using parfor.
|   ├── isdf_kernelg: Generate p_mu(r) and use FFT to get p_mu(g).
|   ├── COmegaCstar: Calculate C * Omega^{-1} * C' with Cauchy integral method,
|       where C_{ij, k} = conj(psi_i(r_k))psi_j(r_k) 
|       and Omega_{ij, ij} = varepsilon_i - varepsilon_j (energies band).
|       CURRENTLY ONLY OMEGA = 0 IS SUPPORTED.
|   └── ellipjc & ellipkkp: a set of ellipses functions to support COmegaCstar, and is directly copied from [1].

├── Utility:
|   ├── testmemory: Test memory usage.

├── Sigma: Contains main function to generate quasiparticle energies
|   ├── gw_cohsex: main function for static approximation.
|   The following are used in GPP calculation.
|   ├── gw_gpp: main function for GPP approximation
|   ├── wpeff: Calculate effective plasmon frequency in GPP calculation.
|   The following are used in full frequency calculation.
|   ├── gw_x: Calculate Hatree--Fock corrolation.
              Currently used in full frequency calculation.
|   ├── gw_fullfreq_cd_res & gw_fullfreq_cd_int: Full frequency calculation of GW using Cauchy integral method.
|   ├── gw_fullfreq_cd_res_isdf & gw_fullfreq_cd_int_isdf: Full frequency calculation of GW using Cauchy integral method and isdf.





[1] https://www.mathworks.com/matlabcentral/fileexchange/1316-schwarz-christoffel-toolbox

