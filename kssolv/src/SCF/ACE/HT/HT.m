function [basis,ev_intp] = HT(psi_g,eigvalues,scfkpts,kpath,mol)
% 1. Hamiltonian transformation to calculate band structure.
%    V(ng,ng,nk)=psi_g(ng,nbnd,nk)*diag(eigvalues(nbnd,nk))*psi_g';
% 2. size(psi_g)=[ng,nbnd,nk]
% 3. size(eigvalues)=[nbnd,nk] or eigvalues=[]
% 4. scfkpts is the kpoints grid of previous scf-calculation, kx varies
%    fastest, kz varies slowest. size(scfkpts)=[nk,3]
% 5. kpath is the kpoint path to be interpolated. size(kpath)=[npath,3]

[~,nbnd,nk]=size(psi_g);
if ~isempty(eigvalues)
    assert(all(size(eigvalues)==[nbnd,nk]))
end
assert(all(size(scfkpts)==[nk,3]))
assert(size(kpath,2)==3)

nk1=mol.nkxyz(1);
nk2=mol.nkxyz(2);
nk3=mol.nkxyz(3);
nr1=mol.n1;
nr2=mol.n2;
nr3=mol.n3;
nr=nr1*nr2*nr3;
nk=nk1*nk2*nk3;

[I,J,K] = ndgrid(0:nr1-1,0:nr2-1,0:nr3-1);
rpts = reshape(cat(4,I,J,K),[],3);
r = rpts./[nr1,nr2,nr3];
scfkpts=scfkpts*mol.supercell'/2/pi;
psi_r=zeros(nr,nbnd,nk);
F=KSFFT(mol);
for ik=1:nk
    psi_r(:,:,ik)=(F'*psi_g(:,:,ik)).*exp(2i*pi*r*scfkpts(ik,:)')*mol.vol/sqrt(nr);
end

psi_r=reshape(psi_r,[nr,nbnd*nk]);

[Q,R,P]=qr(psi_r,0);
threshold=max(abs(diag(R)))*4e-3;
nmu=sum(abs(diag(R))>threshold);
% nmu=48;
invP(P)=1:length(P);
basis=Q(:,1:nmu);
R=R(1:nmu,invP);
Ck=reshape(R,nmu,nbnd,nk);
phi=zeros(nmu,nbnd,nk);
for ik=1:nk
    phi(:,:,ik)=Ck(:,:,ik);
end

Dk=zeros(nmu,nmu,nk);
if isempty(eigvalues)
    for ik=1:nk
        Dk(:,:,ik)=phi(:,:,ik)*phi(:,:,ik)';
    end
else
    gap=4*(max(eigvalues(end,:))-min(eigvalues(end,:)));
    shift=max(eigvalues(end,:));
    power=3;
    ev=trans(gap,power,shift,eigvalues);
    for ik=1:nk
        Dk(:,:,ik)=phi(:,:,ik).*ev(:,ik).'*phi(:,:,ik)';
    end
end

kernel=Finterpolation(Dk,kpath*mol.supercell'/2/pi,mol);
npath=size(kpath,1);

ev_intp=zeros(nbnd,npath);
e_tmp=zeros(nmu,npath);
for i=1:npath
    e=sort(real(eig(kernel(:,:,i))),'ascend');
    ev_intp(1:nbnd,i)=e(1:nbnd);
    e_tmp(:,i)=e;
end
ev_intp = newton_inv(@trans,@dtrans,gap,power,shift,ev_intp);
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

