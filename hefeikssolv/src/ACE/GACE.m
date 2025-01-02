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
[ng,nbnd,nk]=size(zeta_in);

assert(all(size(scfkpts)==[nk,3]))
assert(nk==nkgrid)
assert(size(kpath,2)==3)

[I,J,K] = ndgrid(0:nr1-1,0:nr2-1,0:nr3-1);
rpts = reshape(cat(4,I,J,K),[],3);
r = rpts./[nr1,nr2,nr3];
scfkpts=scfkpts*mol.supercell'/2/pi;
F=KSFFT(mol);

% phase=1 is correct in theory, but the result is worser than phase=0 (which is wrong in theory).
phase=1;
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
ng2=get(F2,'ng');
zeta=zeros(ng2,nbnd,nkgrid);
for ik=1:nk
    zeta(:,:,ik)=F2*((F'*zeta_in(:,:,ik)).*exp(phase*2i*pi*r*scfkpts(ik,:)'));
end

eigval=zeros(nbnd,nkgrid);
for ik=1:nk
    for ib=1:nbnd
        eigval(ib,ik)=norm(zeta(:,ib,ik));
        zeta(:,ib,ik)=zeta(:,ib,ik)/eigval(ib,ik);
        eigval(ib,ik)=eigval(ib,ik)^2;
    end
end
shift=min(eigval,[],'all')*0;
eigval=eigval-shift;

zeta=reshape(zeta,[ng2,nbnd*nk]);
[Q,R,P]=qr(zeta,0);
threshold=max(abs(diag(R)))*1e-5;
n=sum(abs(diag(R))>threshold);
invP(P)=1:length(P);
Q=Q(:,1:n);
R=R(1:n,invP);
Ck=reshape(R,n,nbnd,nkgrid);
Dk=zeros(n,n,nkgrid);
for ik=1:nkgrid
    Dk(:,:,ik)=Ck(:,:,ik).*eigval(:,ik)'*Ck(:,:,ik)';
end
Vk=Finterpolation(Dk,kpath,mol);
zeta=zeros([ng,nbnd*nkgrid,nkpath]);
maxrank=0;
tmp=zeros(nbnd+1,nkpath);
for ik=1:nkpath
    [V,D]=eig(Vk(:,:,ik));
    [d,ind] = sort(diag(D),'descend');
    tmp(:,ik)=d(1:(nbnd+1));
    n = sum(d>max(d)*1e-5);
    maxrank = max(maxrank,n);
    d=d+shift;
    zeta_with_phase=Q*(V(:,ind(1:n)).*sqrt(d(1:n))');
    zeta(:,1:n,ik)=F*((F2'*zeta_with_phase).*exp(phase*-2i*pi*r*kpath(ik,:)'));
end
<<<<<<< HEAD
figure(1)
plot(tmp')
hold on
mol.zeta=zeta(:,1:maxrank,:);
=======
zeta_out=zeta(:,1:maxrank,:);

    
>>>>>>> b5d9dabf0e847d47fcd37d79d9e4105f991a53bf
