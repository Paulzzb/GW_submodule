function Fnl = getFnl(mol,H,X)
% the nonlocal potential contribution to forces
na = sum(mol.natoms);
ecut = mol.ecut;

grid2 = Ggrid(mol,ecut);
gkx   = grid2.gkx;
gky   = grid2.gky;
gkz   = grid2.gkz;

[llsiz,lloffset] = getvnlsize(mol);

if isa(mol,'Crystal')
    nkpts = mol.nkpts;
    Fnl = zeros(na,3);
    for ik = 1:nkpts
        Fnl = Fnl + mol.wks(ik) ...
            *getFnlsingle(H.vnlmatcell{ik},H.vnlsigncell{ik},X{ik});
    end
    if mol.nspin == 2
        for ik = 1:nkpts
            Fnl = Fnl + mol.wks(ik+nkpts) ...
                *getFnlsingle(H.vnlmatcell{ik},H.vnlsigncell{ik},X{ik+nkpts});
        end
    end
else
    Fnl = getFnlsingle(H.vnlmat,H.vnlsign,X);
end
% value should double when nspin == 1
if mol.nspin == 1
    Fnl = Fnl*2;
end

    function Fnlsingle = getFnlsingle(bec,sig,X)
        
        npw = numel(X.idxnz);
        npol = size(X.psi,1)/npw;
        nbeta = sum(llsiz);
        Fnlsingle = zeros(na,3);
        gk = [gkx gky gkz];
      
        if npol == 1                
            becX = bec'*X;
            
            for ipol = 1:3
                dbecX = bec'*(repmat(1i*gk(:,ipol),1,ncols(X)).*X);
                dbecXbecX = real(becX.*conj(dbecX));
                F_tmp = -2*sum(repmat(sig,1,size(becX,2)).*dbecXbecX.*repmat(X.occ',nbeta,1),2);
                for ia = 1:na
                    Fnlsingle(ia,ipol) = sum(F_tmp(lloffset(ia)+(1:llsiz(ia))));
                end
            end
        elseif npol == 2
            % noncolinear case
            X_up = X.psi(1:npw,:);
            X_dw = X.psi(npw+1:end,:);
            
            becX = cell(2,1);
	        becX{1} = bec'*X_up;
	        becX{2} = bec'*X_dw;  
            
            dbecX = cell(2,1);            
            for ipol = 1:3
                gbec = bec.*repmat(-1i.*gk(:,ipol),1,size(bec,2));
                dbecX{1} = gbec'*X_up;
                dbecX{2} = gbec'*X_dw;
                ijs = 0;
                F_tmp = zeros(nbeta,1);
                for is = 1:npol
                    for js = 1:npol
                        ijs = ijs + 1;
                        sig_tmp = diag(sig{ijs});
                        dbecXbecX = conj(dbecX{is}).*becX{js} + conj(becX{is}).*dbecX{js};  
                        F_tmp = F_tmp - sum(repmat(sig_tmp,1,ncols(X)).*dbecXbecX.*repmat(X.occ',nbeta,1),2);                        
                    end
                end
                for ia = 1:na
                    Fnlsingle(ia,ipol) = sum(F_tmp(lloffset(ia)+(1:llsiz(ia))));
                end
                % the additional contribution for full relativistic
                % pseudopotential with multiple projectors
                if mol.lspinorb
                    fac = X.occ;
                    for it = 1:na
                        idx = lloffset(it) + (1:llsiz(it));
                        ijs  = 0;
                        for is = 1:2
                            for js = 1:2
                                ijs = ijs + 1;
                                sig_tmp = sig{ijs}(idx,idx);
                                for n = 1:size(sig_tmp,1)
                                    sig_tmp(n,n) = 0;
                                end
                                
                                sig_up = triu(sig_tmp);
                                sig_dw = tril(sig_tmp);
                                 
                                conj_dbecX = dbecX{is}';
                                conj_becX = becX{is}';
                                Fnlsingle(it,ipol) = Fnlsingle(it,ipol) -fac'*diag(conj_dbecX(:,idx)*sig_up*becX{js}(idx,:) +...
                                    conj_becX(:,idx)*sig_up*dbecX{js}(idx,:));
                                
                                Fnlsingle(it,ipol) = Fnlsingle(it,ipol) - fac'*diag(conj_dbecX(:,idx)*sig_dw*becX{js}(idx,:) +...
                                    conj_becX(:,idx)*sig_dw*dbecX{js}(idx,:));
                            end
                        end
                    end
                end
            end                        
        end                           
    end
end
