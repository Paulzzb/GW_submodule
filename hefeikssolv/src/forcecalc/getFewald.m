function Fewald = getFewald(mol)
% the ewald contribution to the force, very small usually
nel = mol.nel;
alist = mol.alist;
xyzlist = mol.xyzlist;
na = sum(mol.natoms);
ecut2    = mol.ecut2;
% use the numbers of valence electrons in pp
venums = mol.ppvar.venums;
e2 = e2Def();

grid2 = Ggrid(mol,ecut2);
gkk   = grid2.gkk;
gkx   = grid2.gkx;
gky   = grid2.gky;
gkz   = grid2.gkz;

alpha = 1.1;
upperbound = 1;
while upperbound >1e-6
    alpha = alpha-0.1;
    upperbound = e2*nel^2*sqrt(2*alpha/2/pi)*erfc(sqrt(ecut2/2/alpha));
end

if abs(gkk(1)) <= 1e-8
    idxg = 2:length(gkk);
else
    idxg = 1:length(gkk);
end

strf = exp( -1i.*[gkx(idxg) gky(idxg) gkz(idxg)]*xyzlist' );
aux = sum( conj(strf)*diag(arrayfun(@(x)venums(x),alist) ), 2 ) ...
    .*exp(-gkk(idxg)/alpha/4)./gkk(idxg);

Fewald = zeros(na,3);
for it = 1:na
        arg = [gkx(idxg) gky(idxg) gkz(idxg)]*xyzlist(it,:)';
        sumnb = cos(arg).*imag(aux)-sin(arg).*real(aux);
        factor = -venums(alist(it))*2*e2*2*pi/mol.vol;
        Fewald(it,1) = factor*gkx(idxg)'*sumnb;
        Fewald(it,2) = factor*gky(idxg)'*sumnb;
        Fewald(it,3) = factor*gkz(idxg)'*sumnb;
end
if abs(gkk(1)) > 1e-8
    return;
end

rmax = 5/sqrt(alpha);

for iti = 1:na
    for itj = 1:na
        if iti ~= itj
            dxyz = xyzlist(iti,:)-xyzlist(itj,:);
            [r,r2] = rgen(dxyz,rmax);
            if size(r,1) == 0
                continue;
            end
            rr = sqrt(r2);
            factor = venums(alist(iti))*venums(alist(itj))...
                *e2./rr.^2.*( erfc(sqrt(alpha)*rr)./rr ...
                + sqrt(8*alpha/2/pi)*exp(-alpha*rr.^2) );
            Fewald(iti,1) = Fewald(iti,1) - r(:,1)'*factor;
            Fewald(iti,2) = Fewald(iti,2) - r(:,2)'*factor;
            Fewald(iti,3) = Fewald(iti,3) - r(:,3)'*factor;
        end
    end
end

    function [r,r2] = rgen(dxyz,rmax)
        at = mol.supercell;
        bg = inv(mol.supercell);
        n1 = floor(norm(bg(1,:))*rmax)+2;
        n2 = floor(norm(bg(2,:))*rmax)+2;
        n3 = floor(norm(bg(3,:))*rmax)+2;
        
        ds = dxyz*bg;
        ds = ds-round(ds);
        dxyz = ds*at;

        [I,J,K] = ndgrid(-n1:n1,-n2:n2,-n3:n3);
                
         r = [I(:),J(:),K(:)]*at-repmat(dxyz,(2*n1+1)*(2*n2+1)*(2*n3+1),1);
         r2 = sum(r.^2,2);
         idx = r2 <= rmax^2 & r2 > 1e-10;
         r = r(idx,:);
         r2 = r2(idx);
         % maybe no use
         [r2,id] = sort(r2);
         r = r(id,:);  
    end

end
