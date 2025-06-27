function [X,lambada] = davidson_qe(mol,H, X0, tol,maxit,iter,h_diag,s_diag)
% precondition comprises Ekin,Vloc and Vnlc
    %diagt = tic;
    %fprintf("Start davidson diagonalization, convergence threshold is set to %20.13e\n"...
        %,tol);
    npol = 1 + mol.noncolin;
    % set parameters 
    cgoptions = struct;
    cgoptions.tol = tol;
    cgoptions.dav_iter=0;
    cgoptions.maxiter = maxit;
    cgoptions.h_diag = h_diag;
    cgoptions.s_diag = s_diag;
    cgoptions.notcnv = size(X0.psi,2);
    cgoptions.npol = npol;
    cgoptions.btype = ones(1,mol.nbnd);
    
    if ~isa(mol,'Crystal')&&mol.lsda
        e = 2*H.eband((1:mol.nbnd)+mol.nbnd*(X0.ispin-1));
    else
        e = 2*H.eband;
    end
    conv = false;
    ntry = 0;
    while ~conv
        cgoptions.lrot = (iter == 1);
        [X0, e, cgoptions] = diag_david(H,X0,e,cgoptions);
        ntry = ntry + 1;
        if ntry>5 || cgoptions.notcnv == 0
            conv = true;
        end
    end
    % output
    X = X0;
    lambada = e/2;
    %fprintf('avg number of iteration is %d\n', cgoptions.dav_iter);
    %fprintf('Total time for diag = %20.4e\n',toc(diagt));
end
    
    
    
    
