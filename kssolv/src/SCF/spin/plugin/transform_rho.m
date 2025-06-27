function rho = transform_rho(mol,rho,tag)
%transform rho from sum and difference to up and down 
%or inversely
%transform rho from real space to reciprocal space
%or inversely
    rhotmp = rho;
    n123 = mol.n1*mol.n2*mol.n3;
    nspin = mol.nspin;
    if strcmp(tag,'sum2updw')
        rho{1} = (rhotmp{1}+rhotmp{2})/2;
        rho{2} = (rhotmp{1}-rhotmp{2})/2;
    elseif strcmp(tag ,'updw2sum')
        rho{1} = (rhotmp{1}+rhotmp{2});
        rho{2} = (rhotmp{1}-rhotmp{2});
    elseif strcmp(tag,'real2rec')||strcmp(tag , 'rec2real')
        grid = Ggrid(mol,mol.ecut2);
        idxnz2 = grid.idxnz;
        if strcmp(tag,'real2rec')
            if nspin == 2||nspin == 4
                for i = 1:nspin
                    rho{i} = fft3(rho{i})/n123;
                    rho{i} = rho{i}(idxnz2);
                end  
            else
                rho = fft3(rho)/n123;
                rho = rho(idxnz2);
            end
        else
            rhotmp = zeros(mol.n1,mol.n2,mol.n3);
            if nspin == 2||nspin == 4
                for i = 1:nspin
                    rhotmp(idxnz2) = rho{i};
                    rho{i} = ifft3(rhotmp)*n123;
                end
            else
                rhotmp(idxnz2) = rho;
                rho = ifft3(rhotmp)*n123;
            end
        end
    end
end

