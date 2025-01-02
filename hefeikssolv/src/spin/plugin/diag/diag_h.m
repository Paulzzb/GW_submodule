function  [ew,vc] = diag_h(hc,sc,nbase,nvec)

eigsopts.isreal = false;
eigsopts.maxit = 300;
eigsopts.tol = 1e-20;
hc = hc(1:nbase,1:nbase);
sc = sc(1:nbase,1:nbase);
[V,D,flag] = eigs(hc,sc,nvec,'SR',eigsopts);
d = real(diag(D));
[sd,id] = sort(d);
ew = sd;
vc = V(:,id);
if flag~=0
    %fprintf('Convergence not reached in eigs! pause...\');
    pause;
end
end
