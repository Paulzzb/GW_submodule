function kernel = Poisson_solver(mol,zeta,exxgkk)
% Poisson solver for calculating V|zeta(\vr)\rangle
% where V(\vr,\vr^{\prime}) is the Coulomb operator
% The original implementation is kernel=F2'*(exxgkk.*(F2*zeta));
% which is easy but much slower

n1 = mol.n1;n2 = mol.n2;n3 = mol.n3;
n123 = n1*n2*n3;
nu=size(zeta,2);

rhoc=zeros(n1,n2,n3);
vc=zeros(n1,n2,n3);
kernel=zeros(n123,nu);
exxgkk3D=reshape(exxgkk,n1,n2,n3);
for i = 1:nu
    rhoc=reshape(zeta(:,i),n1,n2,n3);
    rhoc=fftn(rhoc);
    vc=rhoc.*exxgkk3D;
    vc=ifftn(vc);
    kernel(:,i)=vc(:);
end

end

