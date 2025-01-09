function [GWinput, options] = qe2GW(options_in)
% function [GWinput, options] = qe2GW(options_in)
  
  % Later need to modify
  if ~isfield(options_in, 'amin')
    amin = 5;
  end
  if ~isfield(options_in, 'inputfile')
    error("Provide output of qe into options.inputfile");
  end

  % READ information from inputfile.
  qepath = options_in.inputfile;
  [sys, options__] = read_qe_gw_bgw(qepath);
  nkpts = sys.nkpts; % Currently, only consider nkpts = 1;



  GWinput = GWinfo();

  % Read Vxc, save it.
  xml=[qepath,'/vxc.dat'];
  fid = dlmread(xml);
  [vxcrow, vxccol] = size(fid);
  step = fid(1,4);
  vxc.kpoints = [];
  vxc.value = [];
  for i =1:step+1:vxcrow-step
  vxc.kpoints=[vxc.kpoints;fid(i,1),fid(i,2),fid(i,3)];
  a=[];
  for j=i+1: 1 :i+step
    a=[a,fid(j,3)];
  end
  vxc.value=[vxc.value; a];
  end
  vxc.value=vxc.value.';
  

  GWinput.Z = options__.X0.wavefuncell{1}.psi;
  GWinput.ev = options__.ev;
  GWinput.Vxc = vxc.value;

  occstatus = options__.X0.wavefuncell{1}.occ;
  GWinput.nv = find((occstatus >= 1-1e-4), 1, 'last');
  GWinput.ne = sys.nel;
  GWinput.vol = sys.vol;
  GWinput.ntot = sys.nbnd;
  GWinput.qk = sys.kpts;
  GWinput.supercell = sys.supercell;
  % Check if need to be squared 
  GWinput.bdot = (2 * pi * inv(sys.supercell)').^2;
  
  % Deal with the grid
  nkpts = length(options__.mill);
  if nkpts == 1
    gvectmp = gvec(); % Construct a empty gvec.
    idxnz_qe = options__.X0.wavefuncell{1}.idxnz;
    mill = options__.mill{1};
    gvectmp.components = mill;
    gvectmp.idxnz = idxnz_qe{1};
    gvectmp.fftgrid = [sys.n1, sys.n2, sys.n3]; 
    gvectmp.nfftgridpts = prod(gvectmp.fftgrid);
    ng = length(idxnz_qe{1});
    gvectmp.ng = ng; 
    GWinput.coulG0 = 8.0 * pi * sys.ecut^2 / 2;
    GWinput.idxnz = gvectmp.idxnz;
    GWinput.gvec = gvectmp;
    GWinput.gvec2;
    GWinput.gvecrho; 
    % Calculate the coulomb potential
    coulG = zeros(ng, 4);
    coulG(:, 1:3) = mill;
    Creci = 2 * pi * inv(GWinput.supercell)';
    kkxyz = double(mill) * Creci;
    gkk = sum(kkxyz.^2, 2);
    for j = 1:ng
      if ( abs(gkk(j)) ~= 0 )
        coulG(j,4) = 8.0*pi/(gkk(j));
        coulG(j,4) = coulG(j,4)*(1-cos(sqrt(gkk(j))*amin));
      else
        %coulG(j) = 4.0*pi*amin^2/2;
        coulG(j,4) = 0.0;
      end
    end 
    GWinput.coulG = coulG;
  else
    ;
  end
 
end % EOF