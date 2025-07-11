function [X,ev] = KSeigs(mol, H, X, prec, options)
%
% usage: [X,ev] = KSeigs(mol, H, X, prec, options);
%
% purupse: update the wavefunctions by computing the invariant
%          subspace associated with the lowest eigenvalues of H
%
eigmethod = options.eigmethod;
verbose   = options.verbose;
if (any(strcmp(options.verbose,{'off';'OFF';'Off'})))
    verbose = 0;
else
    verbose = 1;
end
ncol = ncols(X);

switch lower(eigmethod)
    case {'lobpcg'}
        cgtol = options.cgtol;
        maxcgiter = options.maxcgiter;
        [X, ev, ~, ~] = lobpcg(H, X, prec, cgtol, maxcgiter,verbose);
    case {'davidson'}
        davtol = options.cgtol;
        davsteps = options.davsteps;
        davcycs = options.maxcgiter;
        [X, ev, ~, ~] = davidson(H, X, prec, davtol, davsteps, davcycs, verbose);
    case {'eigs'}
        eigstol = options.eigstol;
        maxeigsiter = options.maxeigsiter;
        [X, ev] = diagbyeigs(H, ncol, eigstol, maxeigsiter);
    case {'chebyfilt'}
        v0 = genX0(mol,1);
        degree = options.degree;
        [T,~,~]=lanczos(H,v0,2*ncol);
        d = sort(real(eig(T)));
        lb = d(ncol+1);
        ub = 1.01*d(2*ncol);
        fprintf('lb = %11.3e, ub = %11.3e, degree = %d\n', ...
            lb, ub, degree);
        Y = chebyfilt(H,X,degree,lb,ub);
        [X,~]=qr(Y,0);
    otherwise
        disp('Unknown method for diagonalizing H! Use eigs');
        [X, ev] = diagbyeigs(H, ncol, eigstol, maxeigsiter);
end
%
%if ( verbose ==1 )
%    HX = H*X;
%    G = X'*HX;
%    R = HX-X*G;
%    ev = sort(real(eig(G)));
%    fprintf('\n');
%    ncol = size(X.psi,2);
%    resnrm = zeros(1,ncol);
%    for j = 1:ncol
%        resnrm(j) = norm(R(:,j))/abs(ev(j));
%        fprintf('eigval(%2d) = %11.3e, resnrm = %11.3e\n', ...
%            j, ev(j), resnrm(j));
%    end
%end
