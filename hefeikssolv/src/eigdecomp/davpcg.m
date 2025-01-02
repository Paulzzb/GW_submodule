function [X,lambda] = davpcg(H, X0, prec, tol, maxit, verbose)

   if (tol < 1e-3)
      [X,lambda,LVEC,RVEC] = lobpcg(H, X0, prec, tol, maxit, verbose);
   else
      [X,lambda] = davidson2(H, X0, prec, tol, maxit, verbose);
   end
end
