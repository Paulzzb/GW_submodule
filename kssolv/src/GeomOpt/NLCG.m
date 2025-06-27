function [mol, FORCE, Etot] = NLCG(mol, ksopts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Amartya S. Banerjee and Chao Yang
% Computational Research Division, Lawrence Berkeley National Lab
%
% Updated by Linhao Wang on 2025/3/22
% Note: several problems still need to be sovled.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

natoms = sum(mol.natoms);

imax = ksopts.maxrelaxiter ;
jmax = 6 ;
n = 30 ;
tol1 = ksopts.relaxtol ;
tol2 = 1e-8 ;
sigma_0 = 0.02 ;

xyzlist = mol.xyzlist;
x = xyzlist(:);

f = @(x)ksefg(x,mol,ksopts) ;

fcalls = 0;

i = 0;
k = 0;

fprintf('\n \n Going for force evaluation type 1 (outside i loop). \n\n');
[Etot,g] = f(x) ; fcalls = fcalls + 1;
fprintf('\n \n Force evaluation type 1, norm(g) = %-16.10e.', norm(g));
fprintf('\n Maximum force component = %-16.10e.\n\n', max(abs(g)));

r = -g ;
M = ones(size(x)) ; % Need to set REAL PRECONDITIONER
s = r./M ;
d = s ;
delta_new = dot(r,d) ;

if imax == 0
    fprintf('\n \n Returning after single shot calculation ... \n\n');
    return;
end

while i<imax && delta_new>tol1
    
    j = 0 ;
    delta_d = dot(d,d) ;
    alph = -sigma_0 ;
    
    fprintf('\n \n Going for force evaluation type 2 (inside i loop, before j loop).');
    fprintf('\n Iteration i = %d.\n\n', i);
    [~,g] = f(x+sigma_0*d); fcalls = fcalls + 1;
    fprintf('\n \n Force evaluation type 2, norm(g) = %-16.10e \n\n', norm(g));
    
    eta_prev = dot(g,d) ;
    
    while (j < jmax) && (alph*alph*delta_d > tol2)
        
        fprintf('\n \n Going in for force evaluation type 3 (inside j-loop).');
        fprintf('\n Iteration i = %d, j = %d.\n\n', i, j);
        
        [~,g] = f(x); fcalls = fcalls + 1;
        fprintf('\n \n Force evaluation type 3, norm(g) = %-16.10e \n\n', norm(g));
        
        eta = dot(g,d) ;
        alph = alph*(eta/(eta_prev-eta)) ;
        x = x + alph*d ;
        eta_prev = eta ;
        j = j + 1 ;
    end
    
    fprintf('\n \n Going in for force evaluation type 4 (inside i loop, after j-loop).');
    fprintf('\n Iteration i = %d.\n\n', i);
    
    [Etot,g] = f(x); fcalls = fcalls + 1;
    fprintf('\n \n Force evaluation type 4, norm(g) = %-16.10e \n\n', norm(g));
    
    r = -g;
    delta_old = delta_new ;
    delta_mid = dot(r,s) ;
    
    M = ones(size(x)) ; % Need to set REAL PRECONDITIONER
    
    s = r./M ;
    delta_new = dot(r,s) ;
    bet = (delta_new-delta_mid)/delta_old ;
    k = k + 1 ;
    
    if (k==n || bet <=0)
        d = s ;
        k = 0 ;
    else
        d = s + bet*d ;
    end
    
    i = i + 1;
    
    fprintf('\n\n Outer CG iteration no. %d completed. \n',i);
    fprintf('\n No. of function (SCF) evaluation calls so far = %d. \n', ...
        fcalls);
    fprintf('\n Current force norm = %-16.10e.\n',norm(g));
    fprintf('\n Maximum force component = %-16.10e.\n\n', max(abs(g)));
    
end

FORCE = reshape(-g,[natoms 3]);

fprintf('\n Total no. of outer NLCG iterations = %d. ', i);
fprintf('\n Total no. of function (SCF) evaluation calls = %d. \n', ...
    fcalls);
fprintf('\n Current force norm = %-16.10e.',norm(g));
fprintf('\n Force vector = \n\n');
disp(FORCE);
fprintf('Maximum force component = %-16.10e. \n\n', max(abs(g)));

xyz =  reshape(x,[natoms 3]);
mol = set(mol,'xyzlist',xyz);

end
