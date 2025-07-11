function [ekin, ord] = sortrx(qq, ng, gvec, sys)

eff_grid = max(gvec(1:ng, :)) - min(gvec(1:ng, :)) + 1;

cholesky = true;
if (cholesky)
    R = chol(sys.bdot);
    qgvec = (qq + gvec) * R.';
    ekin = sum(qgvec.* qgvec,2);
else
    qgvec = (qq + gvec) * 2 * pi * inv(sys.supercell).';
    ekin = sum(qgvec.* qgvec,2);
    % diff nearly zero, two calculations are equivalent.
end

I_G = gvec(:,3) + eff_grid(3) * (gvec(:,2) + eff_grid(2) * gvec(:,1));
rows = 1:size(gvec, 1);
gvec = [gvec ekin I_G rows'];

%Sort by G_square and I_G sequentially
G_rearrange = sortrows(gvec,[4,5]);
G_rearrange(:,[4,5]) = [];

ord = G_rearrange(:, 4);
end