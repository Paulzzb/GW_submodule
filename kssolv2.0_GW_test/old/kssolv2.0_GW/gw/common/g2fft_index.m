function index = g2fft_index(G, gvec)
%G2FFT_INDEX
%   G2FFT_INDEX is used to transform a G vetor in reciprocal lattice
%   to fft vector index. 
G_x = G(:,1);
G_y = G(:,2);
G_z = G(:,3);
% n1 = gvec.n1;
% n2 = gvec.n2;
% n3 = gvec.n3;

index = ((G_x + floor(gvec.n1/2))*gvec.n2 + G_y + floor(gvec.n2/2))*gvec.n3 + G_z + floor(gvec.n3/2) + 1;
end