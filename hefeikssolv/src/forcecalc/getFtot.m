function Ftot = getFtot(mol,H,X,rho)
% the force calculation
% Fnl,Floc,Fscc,Fewald are considerated here
if mol.nspin == 1||mol.nspin == 4
    Fewald = getFewald(mol);
    if ~mol.noncolin
        Floc = getFloc(mol,rho);
    else
        Floc = getFloc(mol,rho{1});
    end
    Fnl = getFnl(mol,H,X);
    Fscc = getFscc(mol,H.dv);
    Ftot = Fewald + Floc + Fnl + Fscc;
elseif mol.nspin == 2
    Fewald = getFewald(mol);
    Floc = getFloc(mol,rho{1}+rho{2});
    if isa(mol,'Crystal')
        Fnl = getFnl(mol,H,X);
    else
        Fnl = getFnl(mol,H,X{1}) + getFnl(mol,H,X{2});
    end
    Fscc = getFscc(mol,H.dv);
    Ftot = Fewald + Floc + Fnl + Fscc;
end
% impose total force = 0
na = length(mol.alist);
total_force = sum(Ftot,1);
Ftot = Ftot - repmat(total_force,na,1)/na;
end