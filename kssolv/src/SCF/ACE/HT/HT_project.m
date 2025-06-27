function [basis,ev_intp] = HT_project(psi_g,eigvalues,scfkpts,kpath,mol)
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
threshold=max(abs(diag(R)))*1e-5;
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
gap=1;
shift=max(eigvalues(end,:));
power=3;
ev=trans(gap,power,shift,eigvalues);
for ik=1:nk
    Dk(:,:,ik)=phi(:,:,ik).*ev(:,ik).'*phi(:,:,ik)';
end

rho=zeros(nmu,nmu,nk);
mu=0.5;
sigma=0.1;
new_ev=0.5*erfc((ev-mu)/sigma);
for ik=1:nk
    rho(:,:,ik)=phi(:,:,ik).*new_ev(:,ik).'*phi(:,:,ik)';
end

Dk=zeros(nbnd,nbnd,nk);
for ik=1:nk
phi_e=phi(:,:,ik).*sqrt(-ev(:,ik).');
[~,~,P]=qr(phi_e',0);
xi=phi_e(P(1:nbnd),:)';
[U,~,V]=svd(xi);
U=U*V';
Dk(:,:,ik)=U'.*ev(:,ik).'*U;
% e1=sort(abs(eig(Dk(:,:,ik))),'ascend');
% e2=sort(abs(ev(:,ik)),'ascend');
% norm(e1(1:nbnd)-e2(1:nbnd))
end

kernel=Finterpolation(Dk,kpath,mol);
npath=size(kpath,1);


ev_intp=zeros(nbnd,npath);
for i=1:npath
    e=sort(real(eig(kernel(:,:,i))),'ascend');
    ev_intp(1:nbnd,i)=e(1:nbnd);
end
ev_intp = newton_inv(@trans,@dtrans,gap,power,shift,ev_intp);
end

function m=distance_matrix(mol)
%Calculate the distance matrix between k-points.
%G=n1*b1+n2*b2+n3*b3, G*ai/2/pi=ni,ni<=|G0|*|ai|/2/pi
B=2*pi*inv(mol.supercell');
kpts=mol.kpts*mol.supercell'/2/pi;
m=zeros(mol.nkpts);
for i=1:mol.nkpts-1
    for j=i+1:mol.nkpts
        G0=mol.kpts(i,:)-mol.kpts(j,:);
        g=kpts(i,:)-kpts(j,:);
        n1_thread=norm(G0)*norm(mol.supercell(1,:))/2/pi;
        n1_min=fix(-n1_thread-g(1));
        n1_max=fix(n1_thread-g(1));
        n2_thread=norm(G0)*norm(mol.supercell(2,:))/2/pi;
        n2_min=fix(-n2_thread-g(2));
        n2_max=fix(n2_thread-g(2));
        n3_thread=norm(G0)*norm(mol.supercell(3,:))/2/pi;
        n3_min=fix(-n3_thread-g(3));
        n3_max=fix(n3_thread-g(3));
        min_dist=norm(g*B);
        for n1=n1_min:n1_max
            for n2=n2_min:n2_max
                for n3=n3_min:n3_max
                    min_dist=min(norm((g+[n1,n2,n3])*B),min_dist);
                end
            end
        end
        m(i,j)=min_dist;
        m(j,i)=min_dist;
    end
end
end

function y=trans(a,n,s,x)
%Transformation function that makes Hamiltonian decays faster
    y=x-s;
    m1=y<=-a;
    m2=y>0;
    y(m1)=y(m1)+a*(1-1/n);
    y(m2)=0;
    y((~m1)&(~m2))=-(-y((~m1)&(~m2))).^n/n/a^(n-1);
end

function y=dtrans(a,n,s,x)
%Derivative of transformation function.
    y=x-s;
    m1=y<=-a;
    m2=y>0;
    y(m1)=1;
    y(m2)=0;
    y((~m1)&(~m2))=(-y((~m1)&(~m2))).^(n-1)*a^(1-n);
end

function x0=newton_inv(f,df,a,n,s,y)
%Newton's method that calculate the inverse of transformation function.
nitermax=100;
dxmax=a;
epsmax=1e-5;
niter=0;
y(y>0)=0;
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
%     disp(niter)
%     disp(max(abs(eps(:))))
end
if max(abs(eps(:)))>epsmax
    fprintf("Warning: newton's method does not converge!\n")
end
end

function data=loadcpl(filename,sz)
    data=fread(fopen(filename),prod(sz)*2,'double');
    data=data(1:2:end)+1i*data(2:2:end);
    data=reshape(data,sz);
end

function data=loadreal(filename,sz)
    data=fread(fopen(filename),prod(sz),'double');
    data=reshape(data,sz);
end

