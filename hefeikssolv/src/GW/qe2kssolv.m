function [chi0, QE_WAV] = qe2kssolv(ksinfo, coulG, omega, eta, nbands)
%
% read wavefunction and Ggrid of Quantum-Espresso
%

%%%Ggrid of kssolv
%C = mol.supercell;
%grid  = Ggrid(mol);
%gkkx  = grid.gkx*C(1,1)/(2*pi);
%gkky  = grid.gky*C(2,2)/(2*pi);
%gkkz  = grid.gkz*C(3,3)/(2*pi);
%gkk   = grid.idxnz;
%
%% construct the discrete Fourier transformation matrix which n123 by ng
%ksinfo.F = KSFFT(mol);
%FF = ksinfo.F;
%[ng,nr]=size(ksinfo.F);
%if (ng ~= length(idxnz))
%   fprintf('ng = %d, length(idxnz) = %d\n', ng, length(idxnz));
%   error('src/GW/gwsetup.m:inconsistent ng and idxnz size');
%end;

%format long

nv    = ksinfo.nv
Z     = ksinfo.Z;
ev    = ksinfo.ev;
F     = ksinfo.F;
vol   = ksinfo.vol;
ntot  = ksinfo.ntot;

im = sqrt(-1);
[ng,nr] = size(F);

idxnz = ng;
%nbands = 6;
fid_1 = fopen ('wfn_qe', 'r');
fid_2 = fopen ('G', 'r');
%fid_3 = fopen ('coulG.dat', 'r')
Format_1 = repmat('%f',1,2);
Format_2 = repmat('%n',1,3);
%Format_3 = repmat('%f',1,4)
QE_WAV_tmp = cell2mat(textscan(fid_1, Format_1, nbands*idxnz));
QE_WAV_tmp = QE_WAV_tmp(:,1) + QE_WAV_tmp (:,2) * im;
QE_WAV = reshape (QE_WAV_tmp,[idxnz,nbands]);
QE_WAV_tmp = QE_WAV;
QE_Ggrid = cell2mat(textscan(fid_2, Format_2, idxnz));
%Ggrid = cell2mat(textscan(fid_3, Format_3, idxnz))
%(m,n) = size (QE_WAV);

% sort wavefunction of Quantum_Espresso with KSSOLV
for i = 1:idxnz
  for j = 1:idxnz
    if coulG(i,1:3) == QE_Ggrid(j,1:3)
      QE_WAV_tmp(i,:) = QE_WAV(j,:);
    end
  end
end
QE_WAV = QE_WAV_tmp;

chi0 = zeros(ng);
fprintf('nr = %d, ng = %d, nbands = %d\n', nr, ng, nbands);

for iv = 1:nbands
  QE_WAV = QE_WAV * sqrt(vol) / (norm(QE_WAV(:,iv)));
end

% Using wavefunctions of QE to calculate M matrix
for iv = 1:nv  % valence band wavefunction
  psiiv = F'*QE_WAV(:,iv);
  for jc = nv+1:nbands  % conduction band wavefunction
    psijc = F'*QE_WAV(:,jc);
    Mr = conj(psiiv).*psijc;
    eden = 1.0/sqrt(2.0*ev(jc)-2.0*ev(iv));
    MR = eden*Mr;
    Mg(:,jc-nv) = F*Mr;
    MG(:,jc-nv) = F*MR;
  end
%  MM = MR;
%  Mg;
  %  save M.dat MM -ascii;
  %fprintf('iv = %d, pairs constructed.\n', iv);
  %fprintf('back to G-space.\n', iv);
  dev = ev(nv+1:nbands) - ev(iv);
  d1 = omega - dev + im*eta;
  d2 = omega + dev - im*eta;
  % 2.0 comes from spin
  D = sparse(diag(1./d1-1./d2));
  % Correct scaling for operators independent of FFT in KSSOLV.
  % Signs consistent with Bruneval (2005)
  %chi0 = chi0 + (conj(Mg)*D)*conj(Mg)'/vol;
  %chi0 = chi0 + (conj(Mg)*D)*conj(Mg)';
  chi0 = chi0 + Mg*(D*Mg') / vol;
  %fprintf('done with inner product.\n', iv);
end;
