function [WFN_FFTgrid] = get_wfn_fftgrid(sys, options)
wfn_box_min = zeros([1 3]);
wfn_box_max = zeros([1 3]);
for ik = 1:sys.nkpts
    [wfn_box_min, wfn_box_max] = get_gvecs_bounds(options.mill{1, ik}, wfn_box_min, wfn_box_max);
end
WFN_FFTgrid = wfn_box_max - wfn_box_min + 1;
end