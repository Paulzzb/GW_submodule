function [h_diag, s_diag] = addnlc(mol,h_diag,sig,vkb)
% add diagonal entries of nonlocal pseudopotential to h_diag
% generate s_diag

na = sum(mol.natoms);
npol = 1 + mol.noncolin;
npw = size(vkb,1);
[llsiz,lloffset] = getvnlsize(mol);

% initialize s_diag to 1
s_diag = ones(npw*npol,1);
ps1 = zeros(2,1);
ps2 = zeros(2,1);
for it = 1:na
    for ih = lloffset(it) + (1:llsiz(it))
        % normal entries
        if mol.noncolin
            ps1(1) = sig{1}(ih,ih);
            ps1(2) = sig{4}(ih,ih);
            ps2(1) = 0;
            ps2(2) = 0;
        else
            ps1(1) = sig(ih);
            ps1(2) = 0;
            ps2(1) = 0;
            ps2(2) = 0;
        end
        vkb_i = vkb(:,ih);
        ar = vkb_i.*conj(vkb_i);
        for ipol = 1:npol                      
            idg = (1:npw)+npw*(ipol-1);
            h_diag(idg) = h_diag(idg) + ps1(ipol)*ar;
%             s_diag(idg) = s_diag(idg) + ps2(ipol)*ar;
        end
        % multiple projection additional entries
        if mol.lspinorb
            for jh = lloffset(it) + (1:llsiz(it))
                if jh~=ih
                    ps1(1) = sig{1}(ih,jh);
                    ps1(2) = sig{4}(ih,jh);
                    ps2(1) = 0;
                    ps2(2) = 0;
                    for ipol = 1:npol
                        vkb_j = vkb(:,jh);
                        ar = vkb_i.*conj(vkb_j);
                        idg = (1:npw)+npw*(ipol-1);
                        h_diag(idg) = h_diag(idg) + ps1(ipol)*ar;
%                         s_diag(idg) = s_diag(idg) + ps2(ipol)*ar;
                    end
                end
            end
        end
    end
end



        
       


    




                    
            
        
    
