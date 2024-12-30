n = 3000;

A = rand(n, n);
eig(A);

AGPU = gpuArray(A);
eig(AGPU);

