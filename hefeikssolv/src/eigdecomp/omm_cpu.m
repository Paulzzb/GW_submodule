function [X,lambda,LVEC,RVEC] = omm_cpu(H, X0, eigub, prec, tol, maxit, verbose);
%
% Usage: [X,lambda,LVEC,RVEC] = omm(H, X0, eigub, prec, tol, maxit, verbose);
%
% Purpose:
%    Compute the invariant subspace associated with the k smallest
%    eigenvalues of the Hamiltonian H using Orbital Minimization, 
%    where k is the number of columns in X0
%
% Input:
%    H     --- a Hamiltonian object
%    X0    --- the initial guess of the wave function in Fourier space
%              (Wavefun object)
%    eigub --- upper bound of the spectrum of H
%    prec  --- preconditioner
%    tol   --- tolarence on the relative residual.
%    maxit --- maximum number of iterations allowed.
%
% Output:
%    X      --- eigenvectors corresponding to the smallest eigenvalues of H
%               (Wavefun object)
%    lambda --- approximations to the smallest eigenvalues of H
%    LVEC   --- eigenvalue history (matrix of size m by ncols, where m is the 
%               total number of iterations and ncols is the number of eigenvalues
%               computed.)
%    RVEC   --- residual norm history (matrix of size m by ncols)
%

%
X      = [];
lambda = [];
LVEC   = [];
%LVEC   = gpuArray(LVEC);
RVEC   = [];
%RVEC   = gpuArray(RVEC);
%
% get size info
%
ncol = ncols(X0); % number of wavefunctions;
if (ncol <=0) 
   fprintf('omm requires at least one initial wave function!\n');
   return;
end
n1=X0.n1;
n2=X0.n2;
n3=X0.n3;
n123 = n1*n2*n3;
%
% orthonormalize the initial wave functions.
%
%X0.psi = gpuArray(X0.psi);
[X,R]=qr(X0,0);
clear X0;
HX = H*X;
%
nconv = 0;
%
iter = 1;
resnrm = ones(ncol,1);  % initialize residual norm
%
% --- M A I N    L O O P ---
%  
while (iter <= maxit & nconv < ncol)
  XTX = X'*X;
  HX = H*X;
  S = X'*HX;
  R = HX-X*S;
  HX = HX - eigub*X;
  S = S - eigub*XTX;
  Grad = 2*HX - X*S - HX*XTX;
  energy = OMMenergy(XTX,S);
  %if (verbose == 1)
  %   fprintf('OMM iter = %3d, energy = %11.3e\n', iter, energy);
  %end 
  %
  % Check for convergence
  %
  for j = 1:ncol 
     resnrm(j) = norm(R(:,j));
     %if (verbose == 1)
     %   fprintf('resnrm(%2d) = %11.3e\n', j, resnrm(j));
     %end
  end
  iconv = find(abs(resnrm)<=tol);
  nconv = length(iconv);
  %
  if (nconv >= ncol) 
     break; 
  end 
  % 
  LVEC(iter,1:ncol) = diag(S)';
  RVEC(iter,1:ncol) = abs(resnrm)';
  %
  % apply the preconditioner prec
  %
  W = -Grad;
  if ( ~isempty(prec) )
     for j = 1:ncol
        W(:,j) = prec.*W(:,j);
     end
  end
  %
  if (iter < 2) 
     %
     % exact line search and update
     %
     [t,X,energy,S,XTX,HX] = linsearchOMM(H,X,W,[],XTX,S,HX,eigub);
  else
     %
     % update the search direction
     %
     beta = real(sum(conj(W).*(W-W0))/sum(conj(W0).*W0));
     W = W+W0*beta;
     %
     % line search along the updated search direction
     %
     [t,X,energy,S,XTX,HX] = linsearchOMM(H,X,W,[],XTX,S,HX,eigub); 
  end
  W0 = W;
  %
  iter = iter + 1;
  %if (verbose == 1)
  %   fprintf('\n'); 
  %end
end
%
% Rayleigh-Ritz for now (can be removed later)
S = X'*HX; S = (S+S')/2;
[Q,D] = eig(S);
[lambda,id] = sort(real(diag(D)));
X = X*Q(:,id);
%X.psi = gather(X.psi);
%lambda = gather(lambda);
%LVEC = gather(LVEC);
%RVEC = gather(RVEC);
end

function [x,Psi,Enew,Hw,Sw,HPsi] = linsearchOMM(H,Psi,D,opt,Sw,Hw,HPsi,shift)
%function [x,Psi,Hw,Sw, Enew] = linMinOMM(H,Psi,D,opt)
% D is the line search direction
% E3 is the energy at the starting point

E3=OMMenergy(Sw,Hw);
Hw2 = D'*(HPsi);
HD = H*D;
Hw3 = D'*(HD-shift*D);

Sw2 = D'*Psi;
Sw3 = D'*D;
c0=E3;
c1=8*trace(Hw2)-4*trace(Sw2*Hw)-4*trace(Sw*Hw2);
c2=4*trace(Hw3)-2*trace(Sw3*Hw)-2*trace(Sw*Hw3)-4*trace(Sw2*Hw2)-4*trace(Sw2'*Hw2);
c3=-4*trace(Sw3*Hw2)-4*trace(Sw2*Hw3);
c4=-2*trace(Sw3*Hw3);
% Important, need to take the real part when the Hamiltonian is a complex
% matrix
c0r = real(c0);
c1r = real(c1);
c2r = real(c2);
c3r = real(c3);
c4r = real(c4);
% solve the quartic function for the global minimizer x
% c0+c1*x+c2*x^2+c3*x^3+c4*x^4
[x fail] = solve_quartic(c0r,c1r,c2r,c3r,c4r);
if (fail)
    error('line search failes');
end
% update
Psi = Psi + D*x;
Hw = Hw + Hw2*x + Hw2'*x + Hw3*x^2;
Sw = Sw + Sw2*x + Sw2'*x + Sw3*x^2;
HPsi = HPsi+HD*x;
Enew = c4*x^4+c3*x^3+c2*x^2+c1*x+c0;

end

function E = OMMenergy(Sw,Hw)
E = trace((2*speye(size(Sw))-Sw)*Hw);
end


function [x_min fail] = solve_quartic(c0,c1,c2,c3,c4)
fail=0;
a=3*c3/(4*c4);
b=2*c2/(4*c4);
if (abs(b)>=1e11 | abs(c4)<=1e-11)
    x_min=-0.5*c1/c2;
else
    d=c1/(4*c4);
    
    Q=(a^2-3*b)/9;
    R=(2*a^3-9*a*b+27*d)/54;
    if (R^2<Q^3)
        theta=acos(R/sqrt(Q^3));
        t=zeros(3,1);  z=t;
        t(1)=-2*sqrt(Q)*cos(theta/3)-a/3;
        t(2)=-2*sqrt(Q)*cos((theta+2*pi)/3)-a/3;
        t(3)=-2*sqrt(Q)*cos((theta-2*pi)/3)-a/3;
        z =c4*t.^4 +c3*t.^3 +c2*t.^2 + c1*t + c0;
        if (c4>0)
            if (min(z(1)>=z(2),z(1)>=z(3)))
                x_order=[1 2 3];
            else if (z(2)>z(3))
                    x_order=[2 3 1];
                else
                    x_order=[3 1 2];
                end
            end
            if (0<=t(x_order(1)) & t(x_order(2))<=t(x_order(1)))
                x_min=t(x_order(2));
            else
                x_min=t(x_order(3));
            end
        else
            if (min(z(1)<=z(2),z(1)<=z(3)))
                x_min=t(1);
            else if (z(2)<z(3))
                    x_min=t(2);
                else
                    x_min=t(3);
                end
            end
        end
    else
        S=-(sign(R)+(R==0))*(abs(R)+sqrt(R^2-Q^3))^(1/3);
        if (S==0)
            U=0;
        else
            U=Q/S;
        end
        x_min=(S+U)-(a/3);
        if (c4<0)
            fail=1;
        end
    end
end
end


