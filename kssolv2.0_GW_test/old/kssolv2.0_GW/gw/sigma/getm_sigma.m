function[gme]=getm_sigma(in, nn, wfnkq, wfnk, iq, ispin, fbz, sys)
%
%% initailization
n1 = sys.n1;
n2 = sys.n2;
n3 = sys.n3;

gme = zeros(fbz.nmtx(iq), sys.nspinor);
%% gvec to fft grid, only translate
g_kq = wfnkq.mill;
g_k = wfnk.mill;
g_q = fbz.mtx{:, iq};
bidx_kq = g2fft_grid(g_kq,n1,n2,n3);
bidx_k = g2fft_grid(g_k,n1,n2,n3);
bidx_q = g2fft_grid(g_q,n1,n2,n3);
%% put into fft_grid
for ispinor = 1 : sys.nspinor
    Nfft1 = zeros(n1,n2,n3);
    Nfft2 = zeros(n1,n2,n3);
    wfnk_in = wfnk.psi{ispin, ispinor}(:,in);
    idx_k = sub2ind(size(Nfft1), bidx_k(:, 1), bidx_k(:, 2), bidx_k(:, 3));
    Nfft1(idx_k) = wfnk_in;
    
    wfnkq_nn = wfnkq.psi{ispin, ispinor}(:,nn);
    idx_kq = sub2ind(size(Nfft2), bidx_kq(:, 1), bidx_kq(:, 2), bidx_kq(:, 3));
    Nfft2(idx_kq) = wfnkq_nn;
    %% inverse fft
    Nfft1 = fftn(Nfft1);
    Nfft2 = fftn(Nfft2);
    %% conjg and multiply
    Nfft1 = conj(Nfft1);
    Nfft2 = Nfft1.*Nfft2;
    %% fft and pick out
    Nfft2 = fftn(Nfft2);
    idx_q = sub2ind(size(Nfft2), bidx_q(:, 1), bidx_q(:, 2), bidx_q(:, 3));
    gme_tmp = Nfft2(idx_q);
    %% scale
    scale = 1/(n1*n2*n3);
    gme(:, ispinor) = gme_tmp * scale;
end
gme = sum(gme, 2);