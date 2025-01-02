function Mk = Finterpolation(M0, kpath, mol)
% Fourier interpolate k-dependent matrix M(:,:,ik) in corse grid to particular k-points path.
% M is 3D array, its last dimension is k-point index.
% The corse grid is mol.scfkpts, which is generated according to mol.nkxyz.
% k-points path is mol.nkpts.

M = reshape(M0,[size(M0,[1,2]),mol.nkxyz]);
nkgrid = prod(mol.nkxyz);
%nkgrid = size(mol.scfkpts,1);
npath = size(kpath,1);

M = fft(M,[],3);
M = fft(M,[],4);
M = fft(M,[],5);
M = reshape(M, [size(M,[1,2]),nkgrid]);

% H0=loadgdb('H0',[36,729,729]);
% H1=loadgdb('H1',[36,729,729]);
% t=0;
% for i =1:36
%     t=max(t,norm(squeeze(H1(i,:,:))-M(:,:,i)/36));
% end

Mk = zeros([size(M,[1,2]),npath]);
n1=mol.nkxyz(1);n2=mol.nkxyz(2);n3=mol.nkxyz(3);
[I,J,K] = ndgrid((0:n1-1)-((0:n1-1) >= n1/2)*n1, ...
    (0:n2-1)-((0:n2-1) >= n2/2)*n2, ...
    (0:n3-1)-((0:n3-1) >= n3/2)*n3);
rpts = reshape(cat(4,I,J,K),[],3);

for ik = 1:npath
    for ir = 1:size(rpts,1)
        Mk(:,:,ik) = Mk(:,:,ik) + exp(2j*pi*rpts(ir,:)*kpath(ik,:)')/nkgrid*M(:,:,ir);
    end
end

% Important!
% Tackle Nyquist frequency when nkx/nky/nkz is even.
for ik = 1:npath
    Mk(:,:,ik)=(Mk(:,:,ik)+Mk(:,:,ik)')/2;
end
