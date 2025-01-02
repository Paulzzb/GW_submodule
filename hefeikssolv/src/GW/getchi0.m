function chi0 = getchi0(ksinfo, omega, eta, nbands);
%
% chi0 = getchi0(ksinfo, omega, eta, nbands);
%
% LL: Revise normalization in an equivalent but consistent form. In
% particular, the normalization of Fourier and real space vectors are
% always preserved. The FFT-specific factors are implemented in-line
% associated with F and F' operators.  (11/28/2015)
%
% NOTE: This implementation does assume that wavefunctions are real, so
% that 
%
% M_{ij}(G) = M_{ji}(G)
%
% The definition of M is slightly different from that in BerkeleyGW but
% is consistent with Bruneval (2005)
%
% ksinfo is a structure that contains ground state calculation results
% and parameters of the molecules
%   nv      --- number of valence states
%   Z       --- contains the eigenvecgtors from the KSDFT calculation
%   ev      --- contains the corresponding eigenvalues
%   vol     --- volume of the unit cell
%   ntot    --- total number of grid points on which the wavefunction is sampled
% 
% omega   --- frequency at which chi0 is being evaluated
% eta     --- Lorentzian broadening factor that turns a Dirac delta into a 
%             a smoother peak
% nbands ---  the total number of bands (KS orbitals) allowed in the chi0 summation. This includes both the occupied and the empty bands

format long

nv   = ksinfo.nv;
Z    = ksinfo.Z;
ev   = ksinfo.ev;
F    = ksinfo.F;
vol  = ksinfo.vol;
ntot = ksinfo.ntot;
%
% 1:nv = occupied, nv+1:n = unoocupied
im = sqrt(-1);
[ng,nr]=size(F);
%%%ZR = F*Z; % transform all wavefunctions to real space 
% need a lot of memory  
% if not enough memory, need to add an extra loop
% to do things on the fly
chi0 = zeros(ng);
if (nargin > 9)
   nbands = ng;
end;
fprintf('nr = %d, ng = %d, nbands = %d\n', nr, ng, nbands);

%nbands = 5;
% Normalize the wavefunction in the Fourier space
% 1/Vol sum_G |Psi(G)|^2 = 1
for iv = 1 : nbands
%  Z(:,iv) = Z(:,iv) / sqrt(vol/ntot) / (norm(Z(:,iv)));
  Z(:,iv) = Z(:,iv) * sqrt(vol) / (norm(Z(:,iv)));
end

%for iv1 = 1:nv
%  psiiv1 = F'*Z(:,iv1);
%  for iv2 = 1:nv
%    psiiv2 = F'*Z(:,iv2);
%    M=F*(conj(psiiv1).*psiiv2)
%  end
%end

for iv = 1:nv
  % Specific scaling due to FFT in KSSOLV.
  psiiv = F'*Z(:,iv);
  for jc = nv+1:nbands
    % Specific scaling due to FFT in KSSOLV.
    psijc = F'*Z(:,jc);
    Mr = conj(psiiv).*psijc;
%    eden = sqrt(1/(2*ev(jc)-2*ev(iv)));
%    MR = eden*Mr;
    % Specific scaling due to FFT in KSSOLV.
    %Mg(:,jc-nv) = conj(F*Mr);
    Mg(:,jc-nv) = F*Mr;
%    MG(:,jc-nv) = F*MR;
  end; 
%  MM = MG;
  Mg;
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
%2.0*ev(1:5);
%MM = MG
%%save M.dat MM -ascii;
%chi0_diag = diag(chi0);
%save chi0.dat chi0_diag -ascii
