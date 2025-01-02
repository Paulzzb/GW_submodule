function ksinfo = gwsetup(mol,amin,nv,nvbands,ncbands,NBANDS)
% usage: set up a bunch of things for GW calculation
% ksinfo = gwsetup(mol)
% ksinfo is a structure that contains the following entries:
%
% nv    number of occupied states (exclude spin)
% Z     contains all eigenvectors of the KS Hamiltonian
% ev    contains the corresponding KS eigenvalues
% F     is the Fourier transform object (non-square) that maps planewave coefficients 
%       to Kohn-Sham orbitals in reals space
% coulG is the Coulomb potential in Fourier space (diagonal of the Coulomb matrix)
% vol   volume of the supercell
% ntot  Total number of grid points in 3D
%
%

format long;

C = mol.supercell;

grid = Ggrid(mol);
gkkx = grid.gkx*C(1,1)/(2*pi);
gkky = grid.gky*C(2,2)/(2*pi);
gkkz = grid.gkz*C(2,2)/(2*pi);
gkk  = grid.gkk; % get a vector |G|^2 within the Ecut limit
idxnz = grid.idxnz; % get the position of gkk within a cube

% construct the discrete Fourier transformation matrix which n123 by ng
ksinfo.F = KSFFT(mol);
FF = ksinfo.F;
[ng,nr]=size(ksinfo.F);
if (ng ~= length(idxnz)) 
   fprintf('ng = %d, length(idxnz) = %d\n', ng, length(idxnz));
   error('src/GW/gwsetup.m:inconsistent ng and idxnz size');
end;

coulG = zeros(ng,4);
for j = 1:length(idxnz)
   coulG(j,1) = gkkx(j);
   coulG(j,2) = gkky(j);
   coulG(j,3) = gkkz(j);
end

nbands = nvbands + ncbands;
im = sqrt(-1);

if nargin>2
  %read wavefunction grid and energy from quantum espresso
  idxnz = ng;
  fid_1 = fopen ('wfn_qe', 'r');
  fid_2 = fopen ('G', 'r');
  fid_3 = fopen ('ev', 'r');
  fid_4 = fopen ('vxc.dat', 'r');
  Format_1 = repmat('%f',1,2);
  Format_2 = repmat('%n',1,3);
  Format_3 = repmat('%f',1);
  Format_4 = repmat('%d %d %f %f',1);
%
%Using eq.(29) in arXiv:1111.4429v3 to calculate the Coulomb-hole operator(Sigma_CH)
% < n k | \Sigma_{CH} (r, r`; 0) | m k > =
% \frac{1}{2} \sum_{q G G`}
% < n k | e^{i (G - G`) \cdot r} | m k >
% [\eps_{G G`}^{-1} (q; 0) - \delta_{G G`}] v (q + G`)
%
  if nargin < 6  
    QE_WAV_tmp = cell2mat(textscan(fid_1, Format_1, nbands*idxnz));
    QE_WAV_tmp = QE_WAV_tmp(:,1) + QE_WAV_tmp (:,2) * im;
    QE_WAV = reshape (QE_WAV_tmp,[idxnz,nbands]);

%
%Using eq.(28) in arXiv:1111.4429v3 to calculate the Coulomb-hole operator(Sigma_CH)
% < n k | \Sigma_{CH} (r, r`; 0) | m k > =
% \frac{1}{2} \sum_{n''} \sum_{q G G`} < n k | e^{i (q + G ) \cdot r} | n'' k >
% < n'' k | e^{-i (q + G' ) \cdot r} | m k >
% [\eps_{G G`}^{-1} (q; 0) - \delta_{G G`}] v (q + G`)
%

  else
    QE_WAV_tmp = cell2mat(textscan(fid_1, Format_1, NBANDS*idxnz));
    QE_WAV_tmp = QE_WAV_tmp(:,1) + QE_WAV_tmp (:,2) * im;
    QE_WAV = reshape (QE_WAV_tmp,[idxnz,NBANDS]); 
  end
  QE_WAV_tmp = QE_WAV;
  QE_Ggrid = cell2mat(textscan(fid_2, Format_2, idxnz));
  ev = cell2mat(textscan(fid_3, Format_3, nbands));
  Vxc = textscan(fid_4,Format_4,nbands);
  ksinfo.vol = det(get(mol,'supercell'));
  ksinfo.ev = ev;
%  vxc{3}
  ksinfo.Vxc = Vxc{3};
  % sort wavefunction of Quantum_Espresso with KSSOLV
  for i = 1:idxnz
    for j = 1:idxnz
      if coulG(i,1:3) == QE_Ggrid(j,1:3)
        QE_WAV_tmp(i,:) = QE_WAV(j,:);
      end
    end
  end
  QE_WAV = QE_WAV_tmp;
  ksinfo.Z = QE_WAV;
  ksinfo.nv = nv;

else
% run SCF to get ground state information
%
  kssolvpptype('pz-hgh','UPF');
  ksinfo.nv = get(mol,'nel')/2; % number of occupied states
  [mol, H, X, info] = scf(mol);
  Hmat = ham2mat(H); % turn it into a matrix
  Hmat = (Hmat+Hmat')/2; % symmetrize to keep eigenvalues real 
  % (causes problem for inversion symmetry?)
  [ksinfo.Z,D] = eig(Hmat);
  [ksinfo.ev,id]=sort(real(diag(D)));
  %ksinfo.ev = ksinfo.ev*2; % convert to Rydberg
  ksinfo.Z = ksinfo.Z(:,id);
  ksinfo.ntot =  get(mol,'n1')*get(mol,'n2')*get(mol,'n3');
  ksinfo.vol = det(get(mol,'supercell'));
%
end


% construct the bare Coulomb in reciprocal space
%coulG = zeros(ng,4);
%scell = mol.supercell;
%amin = 5.0;
%amin = min([norm(scell(:,1)) norm(scell(:,2)) norm(scell(:,3))])/2; % not quite right
for j = 1:ng
   if ( abs(gkk(j)) ~= 0 )
      coulG(j,4) = 8.0*pi/(gkk(j));
      coulG(j,4) = coulG(j,4)*(1-cos(sqrt(gkk(j))*amin));
   else
      %coulG(j) = 4.0*pi*amin^2/2;
      coulG(j,4) = 0.0;
   end;
end;
%coulG;
%save coulG.dat coulG -ascii;
ksinfo.coulG = coulG;
ksinfo.coulG0 = 8.0*pi*amin^2/2;

%if(0)
%%construct the Coulomb interaction with slab truncation
%
%end
%% construct fine K grid to calculate the Coulomb-hole operator sigma_CH
%mol.ecut = 4*mol.ecut;
%grid_fine = Ggrid(mol);
%gkkx_fine = grid_fine.gkx*C(1,1)/(2*pi);
%gkky_fine = grid_fine.gky*C(2,2)/(2*pi);
%gkkz_fine = grid_fine.gkz*C(2,2)/(2*pi);
%%gkk  = grid_fine.gkk;
%idxnz_fine = grid_fine.idxnz;
%G_fine = zeros(length(idxnz_fine),3);
%for i = 1:length(idxnz_fine)
%  G_fine(i,1) = gkkx_fine(i);
%  G_fine(i,2) = gkky_fine(i);
%  G_fine(i,3) = gkkz_fine(i);
%end
%G_fine;
%mol.ecut = mol.ecut/4;
