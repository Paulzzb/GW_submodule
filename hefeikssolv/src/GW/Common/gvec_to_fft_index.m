function idx = gvec_to_fft_index(g, Nfft)
    idx = g + 1;

    if g(1) < 0
        idx(1) = Nfft(1) + idx(1);
    end
    if g(2) < 0
        idx(2) = Nfft(2) + idx(2);
    end
    if g(3) < 0
        idx(3) = Nfft(3) + idx(3);
    end
end
