function [sx, v1x, v2x] = pbexsr_wrap(rho,grho,omega)

nr = length(rho);
sx = zeros(size(rho));
v1x = zeros(size(rho));
v2x = zeros(size(rho));

for ir = 1:nr
    [sx(ir), v1x(ir), v2x(ir)] = pbexsr(rho(ir), grho(ir), omega);
end
end