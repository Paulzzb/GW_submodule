function W = yambo_w_samp_static()
    %YAMBO_W_SAMP_STATIC  w_samp for static COHSEX (NW = 1, mod_FREQUENCIES.F)

    W = struct();
    W.n_freqs = int32(1);
    W.er = single([0, 0]);
    W.ir = single([0, 0]);
    W.damp_reference = single(0);
    W.dr = single([0, 0]);
    W.per_memstps = single(0);
    W.p = complex(single(0), single(0));
    W.samp_type = '2l';
    W.samp_grid = 'lP';
    W.mpa_solver = 'PT';
    gt = 'ra';
    W.grid_type = [gt blanks(16 - numel(gt))];
end
