function W = yambo_w_samp_static()
    %YAMBO_W_SAMP_STATIC  w_samp for static COHSEX (NW = 1, mod_FREQUENCIES.F)

    W = struct();
    W.n_freqs = int32(1);
    W.er = double([0, 0]);
    W.ir = double([0, 0]);
    W.damp_reference = double(0);
    W.dr = double([0, 0]);
    W.per_memstps = double(0);
    W.p = complex(double(0), double(0));
    W.samp_type = '2l';
    W.samp_grid = 'lP';
    W.mpa_solver = 'PT';
    gt = 'ra';
    W.grid_type = [gt blanks(16 - numel(gt))];
end
