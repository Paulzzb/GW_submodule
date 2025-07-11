function gme = getm_sigma(in, nn, wfnkq, wfnk, iq, ispin, fbz, sys)
    % Initialization
    n1 = sys.n1;
    n2 = sys.n2;
    n3 = sys.n3;
    nmtx = fbz.nmtx(iq);
    nspinor = sys.nspinor;
    
    % FFTW optimization
    persistent fftw_planned
    if isempty(fftw_planned)
        fftw('planner', 'measure');
        fftw_planned = true;
    end
    
    % Preallocate memory
    gme = zeros(nmtx, nspinor);
    
    % Convert gvec to FFT grid
    g_kq = wfnkq.mill;
    g_k = wfnk.mill;
    g_q = fbz.mtx{:, iq};
    bidx_kq = g2fft_grid(g_kq, n1, n2, n3);
    bidx_k = g2fft_grid(g_k, n1, n2, n3);
    bidx_q = g2fft_grid(g_q, n1, n2, n3);
    
    % Precompute linear indices
    idx_k = sub2ind([n1, n2, n3], bidx_k(:, 1), bidx_k(:, 2), bidx_k(:, 3));
    idx_q = sub2ind([n1, n2, n3], bidx_q(:, 1), bidx_q(:, 2), bidx_q(:, 3));
    idx_kq = sub2ind([n1, n2, n3], bidx_kq(:, 1), bidx_kq(:, 2), bidx_kq(:, 3));
    
    % Pre-extract psi data
    psi_k_all = cell(nspinor, 1);
    psi_kq_all = cell(nspinor, 1);
    for ispinor = 1:nspinor
        psi_k_all{ispinor} = wfnk.psi{ispin, ispinor}(:, in);
        psi_kq_all{ispinor} = wfnkq.psi{ispin, ispinor}(:, nn);
    end
    
    % Main loop over spinor components
    for ispinor = 1:nspinor
        % Initialize FFT grids
        Nfft1 = zeros(n1, n2, n3);
        Nfft2 = zeros(n1, n2, n3);
        
        % Fill FFT grids
        Nfft1(idx_k) = psi_k_all{ispinor};
        Nfft2(idx_kq) = psi_kq_all{ispinor};
        
        % Optimized FFT calculations
%         Nfft1 = fftn(Nfft1);
%         Nfft2 = fftn(Nfft2);
%         Nfft2 = fftn(conj(Nfft1) .* Nfft2) / (n1 * n2 * n3); 
        Nfft2 = fftn(conj(fftn(Nfft1)) .* fftn(Nfft2)) / (n1 * n2 * n3);
        
        % Extract results
        gme(:, ispinor) = Nfft2(idx_q);
    end
    
    % Sum over spinor components
    gme = sum(gme, 2);
end