function zeta_out = GACE(mol,zeta_in)
% General ACE operator for band calculation.
% Replace the original zeta to new global zeta in mol.

% scfkpts is the kpoints grid of previous scf-calculation, kx varies
% fastest, kz varies slowest. size(scfkpts)=[nkgrid,3]
scfkpts=mol.scfkpts;
% kpath is the kpoint path to be interpolated. size(kpath)=[nkpath,3]
kpath=mol.kpts;
nk1=mol.nkxyz(1);
nk2=mol.nkxyz(2);
nk3=mol.nkxyz(3);
nr1=mol.n1;
nr2=mol.n2;
nr3=mol.n3;
nr=nr1*nr2*nr3;
nkgrid=nk1*nk2*nk3;
nkpath=size(kpath,1);
nspin=mol.nspin;
npol=(nspin==4)+1;
%[ng,nbnd,nk]=size(zeta_in);
if iscell(zeta_in)
    nk = size(zeta_in,1);
    [ng, nbnd] = size(zeta_in{1});
    zeta_tmp = zeta_in;
    zeta_in = zeros(ng, nbnd, nk);
    for ik = 1:nk
        zeta_in(:,:,ik) = zeta_tmp{ik};
    end
else
    [ng,nbnd,nk]=size(zeta_in);
end

assert(all(size(scfkpts)==[nk,3]))
assert(nk==nkgrid)
assert(size(kpath,2)==3)

[I,J,K] = ndgrid(0:nr1-1,0:nr2-1,0:nr3-1);
rpts = reshape(cat(4,I,J,K),[],3);
r = rpts./[nr1,nr2,nr3];
scfkpts=scfkpts*mol.supercell'/2/pi;
F=KSFFT(mol);

% phase=1 is correct in theory, but the result is worser than phase=0 (which is wrong in theory).
phase=0;
if phase==1
    F2=KSFFT(mol,mol.ecut*5);
    fprintf("The warning of KSFFT can be ignored temporary.\n");
    fprintf("This is not wrong, set ecut larger will make F2*(F2'*zeta(:,:,ik))=zeta(:,:,ik)\n");
    fprintf("TODO: The best choice is to set ecut so that max(|G_new|)=max(|G_old|)+1\n");
elseif phase==0
    F2=F;
else
    error('phase not correct')
end
%ng2=get(F2,'ng');
ng2=npol*length(get(F2,'idxnz'));
zeta=zeros(ng2,nbnd,nkgrid);
for ik=1:nk
    for ipol=1:npol
        ig1=(1:ng/npol)+(ipol-1)*ng/npol;
        ig2=(1:ng2/npol)+(ipol-1)*ng2/npol;
        zeta(ig2,:,ik)=F2*((F'*zeta_in(ig1,:,ik)).*exp(phase*2i*pi*r*scfkpts(ik,:)'));
    end
end

eigval=zeros(nbnd,nkgrid);
for ik=1:nk
    for ib=1:nbnd
        eigval(ib,ik)=norm(zeta(:,ib,ik));
        zeta(:,ib,ik)=zeta(:,ib,ik)/eigval(ib,ik);
        eigval(ib,ik)=eigval(ib,ik)^2;
    end
end
%shift=min(eigval,[],'all')*0;
%eigval=eigval-shift;

zeta=reshape(zeta,[ng2,nbnd*nk]);
[Q,R,P]=qr(zeta,0);
threshold=max(abs(diag(R)))*1e-5;
n=sum(abs(diag(R))>threshold);
invP(P)=1:length(P);
Q=Q(:,1:n);
R=R(1:n,invP);
Ck=reshape(R,n,nbnd,nkgrid);
Dk=zeros(n,n,nkgrid);
ev=zeros(nbnd,nkgrid);
maxrank=0;
for ik=1:nkgrid
    Dk(:,:,ik)=Ck(:,:,ik).*eigval(:,ik)'*Ck(:,:,ik)';
    [V,D]=eig(Dk(:,:,ik));
    [d,ind] = sort(real(diag(D)),'descend');
    n = sum(d>1e-13);
    maxrank = max(maxrank,n);
    ev(1:n,ik)=d(1:n);
end
save('Dk.mat','Dk');
gap=4*(max(ev(1,:))-min(ev(1,:)));
shift=max(ev(1,:));
power=3;
%ev=trans(gap,power,shift,ev);
%for ik=1:nkgrid
 %   Dk(:,:,ik)=V(:,ind(1:maxrank)).*ev(1:maxrank,ik)'*V(:,ind(1:maxrank))';
%end
Vk=Finterpolation(Dk,kpath*mol.supercell'/2/pi,mol);
zeta=zeros([ng,nbnd*nkgrid,nkpath]);
maxrank=0;
for ik=1:nkpath
    [V,D]=eig(Vk(:,:,ik));
    %d = newton_inv(@trans,@dtrans,gap,power,shift,diag(D));
    d=diag(D);
    [d,ind] = sort(d,'descend');
    n = sum(d>max(d)*1e-5);
    maxrank = max(maxrank,n);
    zeta_with_phase=Q*(V(:,ind(1:n)).*sqrt(d(1:n))');
    for ipol=1:npol
        ig1=(1:ng/npol)+(ipol-1)*ng/npol;
        ig2=(1:ng2/npol)+(ipol-1)*ng2/npol;        
        zeta(ig1,1:n,ik)=F*((F2'*zeta_with_phase(ig2,1:n)).*exp(phase*-2i*pi*r*kpath(ik,:)'));
    end
end
zeta_out=zeta(:,1:maxrank,:);
end

function y=trans(a,n,s,x)
%Transformation function that makes Hamiltonian decays faster
    y=x-s;
    m1=y<=-a;
    m2=y>0;
    y(m1)=y(m1)+a/2;
    y(m2)=0;
    yt=y((~m1)&(~m2));
    y((~m1)&(~m2))=a*(exp(-n*n/4)-exp(-n*n*(a+2*yt).*(a+2*yt)/4/a/a))/2/n/sqrt(pi)/erf(n/2)+(a+2*yt).*(erf(n/2)-erf(n*(0.5+yt/a)))/4/erf(n/2);
end

function y=dtrans(a,n,s,x)
%Derivative of transformation function.
    y=x-s;
    m1=y<=-a;
    m2=y>0;
    y(m1)=1;
    y(m2)=0;
    y((~m1)&(~m2))=0.5-erf(n*(0.5+y((~m1)&(~m2))/a))/2/erf(n/2);
end

function x0=newton_inv(f,df,a,n,s,y)
%Newton's method that calculate the inverse of transformation function.
nitermax=100;
dxmax=a;
epsmax=1e-5;
niter=0;
x0=y+s;
eps=1;
while (max(abs(eps(:)))>epsmax) && (niter<=nitermax)
    eps=f(a,n,s,x0)-y;
    dx=-eps./df(a,n,s,x0);
    dx(dx>dxmax)=dxmax;
    dx(dx<-dxmax)=-dxmax;
    x0=x0+dx;
    niter=niter+1;
    eps=f(a,n,s,x0)-y;
end
if max(abs(eps(:)))>epsmax
    fprintf("Warning: newton's method does not converge!")
end
end

