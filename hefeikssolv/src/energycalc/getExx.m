function Exx = getExx(X, Vexx, mol)
    nspin = mol.nspin;
    if isa(X,'Wavefun')||~isa(mol,'Crystal')
        % in spin-unrestricted calculations for molecular
        % X is a 2*1 cell array
        if nspin == 1||nspin == 4
       	    VexxPsi = Vexx(X);
            Exx = real(sum(sum(conj(X).*VexxPsi.*X.occ')));
        elseif nspin == 2
            Exx = 0;
             for is = 1:2
                VexxPsi = Vexx(X{is});
                Exx = Exx + real(sum(sum(conj(X{is}).*VexxPsi.*X{is}.occ')));
             end
        end
    elseif isa(X,'BlochWavefun')
        Exx = 0;
	wks = X.wks;
	for i = 1:mol.nkpts*(1+mol.lsda)
            VexxPsi = Vexx(X{i});
            Exx = Exx + wks(i)*real(sum(conj(X{i}).*VexxPsi.*X{i}.occ','all'));
	end
    else
        error('Unknown wavefunction type.')
    end
end
