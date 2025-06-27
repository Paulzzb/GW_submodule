function mag = calmag(mol,rho)
% Calculate the total magnetization for lsda and noncolinear spin
factor = mol.vol/mol.n1/mol.n2/mol.n3;
if mol.lsda
    drho = rho{1} - rho{2};
    mag.absmag = sum(abs(drho),'all')*factor;
    mag.totmag = sum(drho,'all')*factor;
elseif mol.noncolin
    mag.mx = sum(rho{2},'all')*factor;
    mag.my = sum(rho{3},'all')*factor;
    mag.mz = sum(rho{4},'all')*factor;
    mag.absmag = sum(sqrt(rho{2}.^2+rho{3}.^2+rho{4}.^2),'all')*factor;
else
    mag = [];
end
end

