function [sys, options] = gw_setup(sys, options)
%% kpts and grid
options.wfn_fftgrid = get_wfn_fftgrid(sys, options);
options.fftgrid = [sys.n1 sys.n2 sys.n3];
options.kpts = fix(sys.kpts / sys.bvec * 10^15) / 10^15; % reduce accuracy to prevent truncation errors
%% wfn check
for ik = 1:sys.nkpts
    for ispin = 1:sys.nspin
        checknorm(options.X0.wavefuncell{ik, ispin}.psi, ik, ispin);
    end
end
%% fermi energy
rfermi = 0;
efermi_input = 0;
nband = sys.nbnd;
minband = 1;
search_efermi = 1;
update_efermi = 1;
options.efermi = find_efermi(rfermi, efermi_input, nband, minband, search_efermi, update_efermi, sys, options);
%% others
if (sys.nspin == 2)
    options.nv = sys.nel;
else
    options.nv = sys.nel/2;
end
end