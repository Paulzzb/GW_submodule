function Xout = splitX_to_updw(Xin)
% Transform wave function in spin-unrestricted calculations to
% spin-up and spin-down BlochWavefuns
nkpts = Xin.nkpts;
n1 = Xin{1}.n1;n2 = Xin{1}.n2;n3 = Xin{1}.n3;
idxnz = Xin{1}.idxnz;

Psi = cell(nkpts,1);
wks = zeros(nkpts,1);
Xout = cell(2,1);

for is = 1:2
    for ik = 1:nkpts
        id = ik+nkpts*(is-1);
        Psi{ik} = Xin{id}.psi;
        wks(ik) = Xin.wks(id);
    end

    X_tmp = BlochWavefun(Psi,n1,n2,n3,idxnz,wks);    
    for ik = 1:nkpts
        id = ik+nkpts*(is-1);
        X_tmp.wavefuncell{ik}.occ = Xin{id}.occ;
        X_tmp.wavefuncell{ik}.ispin = is;
    end

    Xout{is} = X_tmp;
end

end

