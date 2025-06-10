function [rhoG, gvecrho] = construct_rhoG(data, config)
% construct_rhoG - construct the electrostatic potential from the charge density
% Output:
%     rhoG: density in reciprocal space
%     gvecrho: reciprocal space grid information, @gvec class.

sys = data.sys;
cutoff = config.CUTOFFS.density_cutoff;


n1 = sys.n1;
n2 = sys.n2;
n3 = sys.n3;
n123 = n1*n2*n3;
nfftgrid = [n1, n2, n3];

gvecinput = [];
gvecinput.n1 = n1;
gvecinput.n2 = n2;
gvecinput.n3 = n3;
gvecinput.ecut = cutoff;
gvecinput.supercell = sys.supercell;
gvecrho = gvec(gvecinput);


fftbox = data.rhor;
scal = sys.vol ./ n123;
fftbox = do_FFT(fftbox, nfftgrid, -1);
rhoG = scal * get_from_fftbox(gvecrho.idxnz, fftbox, nfftgrid);

end % EOF
