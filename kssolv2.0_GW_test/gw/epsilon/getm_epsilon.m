function gme = getm_epsilon(ivin, icin, wfnkq, wfnk, iq, ispin, pol, sys, options)
    %% FFTW optimization
    persistent fftw_planned
    if isempty(fftw_planned)
        fftw('planner', 'measure');
        fftw_planned = true;
    end
    
    %% set up fft_grid
    Nfac = 3;
    for i = 1:3
        Nfft(i) = pol.fftgrid{:, iq}(1, i);
        while (not(check_FFT_size(Nfft(i), Nfac)))
            Nfft(i) = Nfft(i) + 1;
        end
    end

    gme = zeros(pol.nmtx(iq), sys.nspinor);
    %% gvec to fft grid, only translate
    g_kq = wfnkq.mill;
    g_k = wfnk.mill;
    g_q = pol.mtx{:, iq};
    bidx_kq = g2fft_grid(g_kq, Nfft(1), Nfft(2), Nfft(3));
    bidx_k = g2fft_grid(g_k, Nfft(1), Nfft(2), Nfft(3));
    bidx_q = g2fft_grid(g_q, Nfft(1), Nfft(2), Nfft(3));
    
    %% Precompute indices and wavefunction data
    idx_kq = sub2ind([Nfft(1), Nfft(2), Nfft(3)], bidx_kq(:,1), bidx_kq(:,2), bidx_kq(:,3));
    idx_k = sub2ind([Nfft(1), Nfft(2), Nfft(3)], bidx_k(:,1), bidx_k(:,2), bidx_k(:,3));
    idx_q = sub2ind([Nfft(1), Nfft(2), Nfft(3)], bidx_q(:,1), bidx_q(:,2), bidx_q(:,3));
    
    vwfn_all = cell(sys.nspinor, 1);
    cwfn_all = cell(sys.nspinor, 1);
    for ispinor = 1:sys.nspinor
        vwfn_all{ispinor} = wfnkq.psi{ispin, ispinor}(:, ivin);
        cwfn_all{ispinor} = wfnk.psi{ispin, ispinor}(:, icin);
    end    
    %% Main computation loop
    for ispinor = 1:sys.nspinor
        Nfft1 = zeros(Nfft(1), Nfft(2), Nfft(3));
        Nfft2 = zeros(Nfft(1), Nfft(2), Nfft(3));
        
        % Fill FFT grids
        Nfft1(idx_kq) = vwfn_all{ispinor};
        Nfft2(idx_k) = cwfn_all{ispinor};
        
        % FFT computations
%         Nfft1 = fftn(Nfft1);
%         Nfft2 = fftn(Nfft2);
%         Nfft2 = fftn(conj(Nfft1) .* Nfft2) / (Nfft(1)*Nfft(2)*Nfft(3));
        Nfft2 = fftn(conj(fftn(Nfft1)) .* fftn(Nfft2)) / (Nfft(1)*Nfft(2)*Nfft(3));
        
        % Extract and scale results
        gme(:, ispinor) = Nfft2(idx_q);
    end
    
    %% Sum over spinor components
    gme = sum(gme, 2);
    
    %% get eden and multiply
    eval = options.ev(ivin, wfnkq.ikq, ispin);
    econd = options.ev(icin, wfnk.ikq, ispin);
    eden = 1/sqrt(eval-econd);
    gme = gme * eden;
end