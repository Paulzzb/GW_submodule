function[gme] = getm_epsilon(ivin, icin, wfnkq, wfnk, iq, ispin, pol, sys, options)
%% set up fft_grid
Nfac = 3;
for i = 1 : 3
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
%% put into fft_grid
for ispinor = 1 : sys.nspinor
    Nfft1 = zeros(Nfft(1),Nfft(2),Nfft(3));
    Nfft2 = zeros(Nfft(1),Nfft(2),Nfft(3));
    vwfn = wfnkq.psi{ispin, ispinor}(:, ivin);
    idx_kq = sub2ind(size(Nfft1), bidx_kq(:, 1), bidx_kq(:, 2), bidx_kq(:, 3));
    Nfft1(idx_kq) = vwfn;
    
    cwfn = wfnk.psi{ispin, ispinor}(:, icin);
    idx_k = sub2ind(size(Nfft2), bidx_k(:, 1), bidx_k(:, 2), bidx_k(:, 3));
    Nfft2(idx_k) = cwfn;
    %% inverse fft
    Nfft1 = fftn(Nfft1);
    Nfft2 = fftn(Nfft2);
    %% conjg and multiply
    Nfft1 = conj(Nfft1);
    Nfft2 = Nfft1.* Nfft2;
    %% fft and pick out
    Nfft2 = fftn(Nfft2);
    idx_q = sub2ind(size(Nfft2), bidx_q(:, 1), bidx_q(:, 2), bidx_q(:, 3));
    gme_tmp = Nfft2(idx_q);
    %% scale
    scale = 1 / (Nfft(1) * Nfft(2) * Nfft(3));
    gme(:, ispinor) = gme_tmp * scale;
end
gme = sum(gme, 2);
%% get eden and multiply
eval = options.ev(ivin, wfnkq.ikq, ispin);
econd = options.ev(icin, wfnk.ikq, ispin);
eden = 1/sqrt(eval-econd);
gme = gme * eden;

