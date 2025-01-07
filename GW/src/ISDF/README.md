# Structure in source code
```


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
|   ├── ellipjc & ellipkkp: a set of ellipses functions to support COmegaCstar, and is directly copied from [1].
|   └── landen & ellipk: two functions from siglib,
                         used in COmegaCstar.

``````


[1] https://www.mathworks.com/matlabcentral/fileexchange/1316-schwarz-christoffel-toolbox

