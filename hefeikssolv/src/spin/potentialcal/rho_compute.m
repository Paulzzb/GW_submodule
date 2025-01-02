function [rhout,segni] = rho_compute(rho,lsign,ux)
% rotate rho through a unitary matrix to up and down format
% such that we can use LSDA to calculate functional
rhout = cell(2,1);
if lsign
    segni = sign(rho{2}*ux(1)+rho{3}*ux(2)+rho{4}*ux(3));
    amag = sqrt(rho{2}.^2+rho{3}.^2+rho{4}.^2);
    rhout{1} = 0.5*(rho{1}+segni.*amag);
    rhout{2} = 0.5*(rho{1}-segni.*amag);
else
    segni = ones(size(rho{1}));
    amag = sqrt(rho{2}.^2+rho{3}.^2+rho{4}.^2);
    rhout{1} = 0.5*(rho{1}+amag);
    rhout{2} = 0.5*(rho{1}-amag);
end
end