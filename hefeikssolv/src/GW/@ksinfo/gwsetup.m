function ksinfo = gwsetup(ksinfo, mol, options)

% This function sets up a bunch of things for GW calculation
%
% Input:
% mol: with class molecule.
% options.:
%   amin: use to take sphere cutoff, with Ry as its unit.
%   nvbands: number of valence bands to calculate, counted from Fermi level.
%   ncbands: number of conductive bands to calculate, counted from Fermi level.
%   NBANDS: nvbands + ncbands.
%   BSE_options.: not appliable yet.
%   input: 'qe' or 'kssolv'.
%
% Output: 
% ksinfo: structure, some important entries as followed:
%   nv      Number of occupied states (exclude spin)
%   Z       KS--DFT wavefunctions on fourier grid.
%   ev      KS--DFT energies.
%   F       Fourier transform object (non-square) that maps planewave coefficients 
%           to Kohn-Sham orbitals in reals space.
%   coulG   Coulomb potential in Fourier space (diagonal of the Coulomb matrix).
%   coulG0  Coulomb potential at q+G --> 0.
%   vol     volume of the supercell.
%   ntot    Total number of grid points in 3D.
%   Vxc     Exchange and corrolation energies.
%   rpoints   maps between 1-D real grid with 3-D real grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note:
% All calculations is based on Ry. 
% Energies in input file should always transpose to Ry.

% Now only allowed options.input = 'kssolv'. 
% if options.input = 'kssolv'
%   use scf.m to calculate datas of ground states.
%   transform energies from Ha to Ry.

% nameConstants = fieldnames(options.Constant);
% for i = 1:numel(nameConstants)
%     fieldname = nameConstants{i};
%     value = options.Constant.(fieldname);    
%     % Now, you can use the variable "fieldname" and "value" in your function
%     strEval = sprintf('%s = %f;', fieldname, value);
%     eval(strEval);
%     fprintf('Field: %s, Value: %s\n', fieldname, num2str(value));
% end


% Some initialization.
if ~isfield(options, 'options_kssolv')
	default_scf_options = true;
else
	default_scf_options = false;
end


if ~isfield(options, 'coulomb_truncation')
	options.coulomb_truncation = 2;
end



% Constant
im    = sqrt(-1);
ha2ry = 2;
ry2ev = 13.60569253;


% Inputs from options
nv = options.nv;
nc = options.nc;
nbands = nv + nc;
amin = options.amin
ng = 0;
nr = 0;

% Grid maping relation calculation
grid = Ggrid(mol);
C = mol.supercell;
gkkx = grid.gkx * C(1,1) / (2*pi);
gkky = grid.gky * C(2,2) / (2*pi);
gkkz = grid.gkz * C(3,3) / (2*pi);
gkk  = grid.gkk; % get a vector |G|^2 within the Ecut limit
idxnz = grid.idxnz; % get the position of gkk within a cube

% construct F --> @KSFFT 
F = KSFFT(mol);
FFT_grid = [get(mol,'n1'),get(mol,'n2'),get(mol,'n3')];
[ng,nr] = size(F);

if (ng ~= length(idxnz)) 
   fprintf('ng = %d, length(idxnz) = %d\n', ng, length(idxnz));
   error('src/GW/gwsetup.m:inconsistent ng and idxnz size');
end

coulG = zeros(ng, 4);
for j = 1:ng
  coulG(j,1) = int32(gkkx(j));
  coulG(j,2) = int32(gkky(j));
  coulG(j,3) = int32(gkkz(j));
end

%spherical_truncation, ry
switch options.coulomb_truncation
	case 2 % spherical_truncation
    for j = 1:ng
      if ( abs(gkk(j)) ~= 0 )
          coulG(j,4) = 8.0*pi/(gkk(j));
          coulG(j,4) = coulG(j,4)*(1-cos(sqrt(gkk(j))*amin));
      else
          %coulG(j) = 4.0*pi*amin^2/2;
          coulG(j,4) = 0.0;
      end
    end
  case 6 %cell_slab_truncation
    zc = mol.supercell(3,3) / 2;
    %get(mol,'supercell')
    for j = 1:ng
      if ( abs(gkk(j)) ~= 0 )
          coulG(j,4) = 8.0*pi/(gkk(j));
          gkkxy = sqrt(grid.gkx(j)^2+grid.gky(j)^2);
    %      grid.gkz
          coulG(j,4) = coulG(j,4)*(1-exp(-gkkxy*zc)*cos(grid.gkz(j)*zc));
      else
          %coulG(j) = 4.0*pi*amin^2/2;
          coulG(j,4) = 0.0;
      end;
    end;
  otherwise
	  fprintf('options.coulomb_truncation = %d is not supported.\n', ...
		         options.coulomb_truncation);
    error();
end;


ksinfo.F = F;
ksinfo.mol = mol;
ksinfo.coulG = coulG;
ksinfo.coulG0 = 8.0*pi*amin^2/2;
ksinfo.bdot = (2 * pi * inv(ksinfo.mol.supercell)').^2;    
ksinfo.qk = [0.0, 0.0, 0.0];
ksinfo.ntot =  get(mol,'n1') * get(mol,'n2') * get(mol,'n3');
ksinfo.vol = det(mol.supercell);
ksinfo.ne = get(mol, 'nel');
ksinfo.nv = nv; 


switch lower(options.input)
	case 'kssolv'
  % run SCF to get ground state information
  kssolvpptype('pz-hgh', 'UPF');
  if default_scf_options
	  options.options_kssolv = setksopt();
    options.options_kssolv.nc = nc; % get an extra band.
    options.options_kssolv.useace = 1;
    %options._kssolv.what2mix = 'rho';
    %options._kssolv.mixtype = 'broyden';
    options.options_kssolv.eigmethod = 'eigs';
    options.options_kssolv.verbose = 0;
    options.options_kssolv.maxphiiter = 30;
    options.options_kssolv.phitol = 1e-8;
  end
  if isfield(options, 'inputfile')
    load(options.inputfile);
    [~, vxc] = getVhxc(mol, H.rho);
  else
    [mol,H,X0,info] = scf(mol, options.options_kssolv);
    vxc = H.vxc;
  end

% calculate exchange-correlation potential Vxc = <nk|vxc|nk>
  vxc_tmp = reshape(vxc, nr, 1);
  psir = F' * X0.psi;
  Vxc = zeros(nbands, 1);
  for it=1:nbands
    Vxc(it) = sumel(vxc_tmp .* ((abs(psir(:,it))).^2)) * (ksinfo.vol)^2/nr;
  end
  ksinfo.Vxc = Vxc * ha2ry;
  
	if options.frequency_dependence == 1
    rhor = H.rho; 
    atomlist = mol.atoms(mol.alist);
    mol_rho = Molecule('supercell',mol.supercell, 'atomlist', atomlist,...
		                   'xyzlist' ,mol.xyzlist, 'ecut', mol.ecut * 4, ...
											 'n1', mol.n1, 'n2', mol.n2, 'n3', mol.n3, 'nbnd', mol.nbnd);
    grid_rho = Ggrid(mol_rho);
    idxnz_rho = grid_rho.idxnz; 
		Nfft = [mol.n1, mol.n2, mol.n3];
		coulG_rho = [grid_rho.gkx, grid_rho.gky, grid_rho.gkz] * mol_rho.supercell / (2*pi);
		coulG_rho = int64(coulG_rho);
		gindex = 1:length(grid_rho.gkx);
		fftbox2 = rhor;
		fftbox2 = do_FFT(fftbox2, Nfft, -1);
		scale = 1 ./ nr;
    % rhog = ksinfo.vol * scale * get_from_fftbox(length(gindex), coulG_rho(:, 1:3), gindex, fftbox2, Nfft);
    rhog = ksinfo.vol * scale * get_from_fftbox(idxnz_rho, fftbox2, Nfft);
		ksinfo.rho = rhog;
    ksinfo.gvecrho = gvec(mol_rho);
  end
	
%	gvec = struct();
%  if (options.frequency_dependence == 1)
%	  gvec.components = coulG_rho;
%	  gvec.index_vec(grid_rho.idxnz) = (1:grid_rho.ng)';
%	  gvec.ng = grid_rho.ng;
%	else
%		gvec.components = coulG(:, 1:3);
%		gvec.index_vec(grid.idxnz) = (1:grid.ng)';
%		gvec.ng = grid.ng;
%	end
%	gvec.nfftgridpts = nr;
%	gvec.fftgrid = [mol.n1, mol.n2, mol.n3];
%  ksinfo.gvec = gvec;
	ksinfo.gvec = gvec(mol, 'ecut', mol.ecut);
	ksinfo.gvec2 = gvec(mol, 'ecut', mol.ecut2);
  ksinfo.idxnz = idxnz; 
 
  % Other important information of DFT-KSSOLV
  ksinfo.ev = info.Eigvals * ha2ry; % transform Ha to Ry
  ksinfo.Z = X0.psi(:, 1:nv+nc);
% plasma_omega = 4.0 * sqrt(2*pi * ksinfo.nv/ksinfo.vol);
% ksinfo.plasma_omega = plasma_omega;

   
	  
	otherwise
	  error('The case is not supported yet!')
end

% Gives outputs
fprintf('system = %s\n', mol.name);
fprintf('Cutoff of Wavefunctions and dielectric matrix = %d (Ry)\n', 2*mol.ecut);
fprintf('Number of valence bands = %d\n', nv);
fprintf('Number of conduction bands = %d\n', nc);
fprintf('nr = %d, ng = %d, nbands = %d\n', nr, ng, nbands);
fprintf('setting ksinfo done!\n')
fprintf('\n')

return % main function
end % main function






