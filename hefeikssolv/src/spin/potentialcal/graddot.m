function da = graddot(mol,a)
% Calculation of laplace operation of the function whose first 
% derivation is a
n1 = mol.n1;
n2 = mol.n2;
n3 = mol.n3;
ecut2 = mol.ecut2;

grid2 = Ggrid(mol,ecut2);
idxnz2 = grid2.idxnz;
gk2 = cell(3,1);
gk2{1} = grid2.gkx;
gk2{2} = grid2.gky;
gk2{3} = grid2.gkz;
gaig = zeros(n1,n2,n3);
for i = 1 : 3
    ai = complex(a{i});
    aig = fft3(ai);
    aig   = aig(idxnz2);    
    gaig(idxnz2) = gaig(idxnz2)+1i*aig.*gk2{i};       
end
    da = real(ifft3(reshape(gaig,n1,n2,n3)));
end

