function [basis,ev_intp] = WPO(psi_g,eigvalues,scfkpts,kpath,mol)
% 1. Weighted projection operator V
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

% wfc=loadcpl('wfc',size(psi_r));
% tmp=loadgdb('tmp',[729,28,36]);
% t1=zeros(nr,nbnd,nk);
% t1(mol.qeidxnz,:,:)=psi_g;
n_discard=0;
nb_part=size(psi_r,2)-n_discard;
psi_r_part=psi_r(:,n_discard+1:end,:);

%psi_r=reshape(psi_r,[nr,nbnd*nk]);
psi_r_part=reshape(psi_r_part,[nr,nb_part*nk]);

% wfc=loadcpl('wfc',size(psi_r_part));
% rk=min(size(psi_r));
% Q=loadcpl('Q',[size(psi_r,1),size(psi_r,2)]);
% Q=Q(:,1:rk);
% R=loadcpl('R',[rk,size(psi_r,2)]);
% nmu=size(R,1);
% basis=Q(:,1:nmu);

[Q,R,P]=qr(psi_r_part,0);
threshold=max(abs(diag(R)))*1e-5;
nmu=sum(abs(diag(R))>threshold);
invP(P)=1:length(P);
basis=Q(:,1:nmu);
R=R(1:nmu,invP);
Ck=reshape(R,nmu,nb_part,nk);
phi=zeros(nmu,nbnd,nk);
for ik=1:nk
    phi(:,1:n_discard,ik)=basis'*psi_r(:,1:n_discard,nk)*0;
    phi(:,n_discard+1:nbnd,ik)=Ck(:,:,ik);
end

Dk=zeros(nmu,nmu,nk);

if isempty(eigvalues)
    for ik=1:nk
        Dk(:,:,ik)=phi(:,:,ik)*phi(:,:,ik)';
    end
else
    ev=eigvalues;
    sigma=(max(ev(end,:))-min(ev(end,:)));
    mu=min(ev(end,:))+10;
    s=max(ev(end,:));
    ev=(ev-s).*erfc((ev-mu)/sigma)/2;
    for ik=1:nk
        Dk(:,:,ik)=phi(:,:,ik).*ev(:,ik).'*phi(:,:,ik)';
    end
end

% H=loadcpl('H',[nk,rk,rk]);
% for i =1:nk
%     norm(squeeze(H(i,:,:))-Dk(:,:,i))
% end

kernel=Finterpolation(Dk,kpath,mol);
npath=size(kpath,1);

% M=load('M0',[npath,size(psi_r,1),size(psi_r,1)]);
% for i =1:npath
%     norm(squeeze(M(i,:,:))-kernel(:,:,i))
% end

ev_intp=zeros(nbnd,npath);
for i=1:npath
    e=sort(real(eig(kernel(:,:,i))),'ascend');
    ev_intp(1:nbnd,i)=e(1:nbnd);
end
ev_intp = inv_eigvalues(ev_intp);

    function x=inv_eigvalues(y)
        % solve (x-s)erfc((x-mu)/sigma)/2=y
        x=y+s;
        % This inverse is not correct, since the transformation is 
        % not monotonic function if max(y)>s.
        for iter=1:20
            f=(x-s).*erfc((x-mu)/sigma)/2-y;
            %f=f.*(y<0);
            if max(abs(f(:)))<max(abs(x(:)))*1e-6
                break
            end
            df=-(x-s)/sigma/sqrt(pi).*exp(-((x-mu)/sigma).^2)+erfc((x-mu)/sigma)/2;
            dx=-f./df;
            dx(dx>sigma/2)=sigma/2;
            dx(dx<-sigma/2)=-sigma/2;
            x=x+dx;
        end
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
