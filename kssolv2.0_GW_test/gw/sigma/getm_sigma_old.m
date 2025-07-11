function gme = getm_sigma(in, nn, wfnkq, wfnk, iq, ispin, fbz, sys)
    % Initialization
    n1 = sys.n1;
    n2 = sys.n2;
    n3 = sys.n3;
    nmtx = fbz.nmtx(iq);
    nspinor = sys.nspinor;
    
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
    
    % Main loop over spinor components
    for ispinor = 1 : nspinor
        % Initialize FFT grids
        Nfft1 = zeros(n1, n2, n3);
        Nfft2 = zeros(n1, n2, n3);
        
        % Fill FFT grids
        wfnk_in = wfnk.psi{ispin, ispinor}(:, in);
        Nfft1(idx_k) = wfnk_in;
        
        wfnkq_nn = wfnkq.psi{ispin, ispinor}(:, nn);
        idx_kq = sub2ind([n1, n2, n3], bidx_kq(:, 1), bidx_kq(:, 2), bidx_kq(:, 3));
        Nfft2(idx_kq) = wfnkq_nn;
        
        % Perform FFT and multiplication
        Nfft1 = fftn(Nfft1);
        Nfft2 = fftn(Nfft2);
        Nfft2 = conj(Nfft1) .* Nfft2;
        Nfft2 = fftn(Nfft2);
        
        % Extract results and scale
        gme(:, ispinor) = Nfft2(idx_q) / (n1 * n2 * n3);
    end
    
    % Sum over spinor components
    gme = sum(gme, 2);
end