function index = index_mix(j,nspin)
% the index of 2D array for different nspin
index = (1:nspin) + nspin*(j-1);
end