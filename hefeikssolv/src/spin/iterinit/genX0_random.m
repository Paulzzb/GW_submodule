function [X0, H] = genX0_random(mol,H)
% The random generation method of wavefunctions
if isa( mol, 'Crystal' )
    nXcols = mol.nbnd*ones(mol.nkpts,1);
    nkpts = mol.nkpts;
else
    nXcols = mol.nbnd;
    nkpts = 1;
end
n1   = mol.n1;
n2   = mol.n2;
n3   = mol.n3;
nspin = mol.nspin;
grid = Ggrid(mol);
idxnz = grid.idxnz;

if nspin == 2
    Qcell = cell(nkpts*2,1);
else
    Qcell = cell(nkpts,1);
end
tpiba = 2*pi/mol.alat;
info.first = true;
[gkk,id] = sort(grid.gkk);
map = zeros(numel(id),1);
for i = 1:numel(id)
    for j = 1:numel(id)
        if id(j) == i
            map(i) = j;
        end
    end
end
npol = 1+mol.noncolin;
npw = size(idxnz,1);
for ik = 1:nkpts  
  nbnd = nXcols(ik);
  psif = zeros(npw,npol,nbnd);  
  for j = 1:nbnd
      if isa(mol,'Crystal')
          gg = sumel(mol.kpts(ik,:).^2);
      else
          gg = 0;
      end
      for ipol = 1:npol
          for ig = 1:npw
              [rr,info] = randy(info);
              [arg,info] = randy(info);
              arg = 2*pi*arg;
              psif(ig,ipol,j) = complex(rr*cos(arg),rr*sin(arg))/(1+gkk(ig)/...
                  tpiba^2+gg/tpiba^2);
          end
      end
  end

  if ~mol.noncolin
      Qcell{ik} = squeeze(psif(map,1,:));
  else
      psif = psif(map,:,:);
      Qcell{ik} = zeros(npw*npol,nbnd);
      Qcell{ik}(1:npw,1:nbnd) = psif(:,1,:);
      Qcell{ik}(npw+1:end,1:nbnd) = psif(:,2,:);
  end
end

if nspin == 2
     for ik = 1:nkpts  
      psif = zeros(size(idxnz,1),nXcols(ik));  
      for j = 1:nXcols(ik)  
          if isa(mol,'Crystal')
              gg = sumel(mol.kpts(ik,:).^2);
          else
              gg = 0;
          end
          for ipol = 1:npol
              for ig = 1:size(idxnz,1)
                  [rr,info] = randy(info);
                  [arg,info] = randy(info);
                  arg = 2*pi*arg;
                  psif(ig,j) = complex(rr*cos(arg),rr*sin(arg))/(1+gkk(ig)/...
                      tpiba^2+gg/tpiba^2);
              end
          end
      end
      Qcell{ik+nkpts} = psif(map,:);
     end
end
% Rotate wavefunctions
if nspin == 2
    if isa( mol, 'Crystal' )
        X0 = BlochWavefun(Qcell,n1,n2,n3,idxnz,mol.wks,mol.nspin);
        H.eband = cell(nkpts*2,1);
        for ik = 1:nkpts
            [X0{ik},H.eband{ik}] = rotate_wfc(mol, H{ik}, X0{ik});
            [X0{ik+nkpts}, H.eband{ik+nkpts}] = rotate_wfc(mol, H{ik}, X0{ik+nkpts});
        end
    else
        X0 = cell(2,1);
        X0{1} = Wavefun(Qcell{1},n1,n2,n3,idxnz,1);
        X0{2} = Wavefun(Qcell{2},n1,n2,n3,idxnz,2);
        H.eband = zeros(nXcols*2,1);
        [X0{1}, H.eband(1:nXcols)] = rotate_wfc(mol, H, X0{1});
        [X0{2}, H.eband(nXcols+1:end)] = rotate_wfc(mol, H, X0{2});
    end
else
    if isa( mol, 'Crystal' )
        X0 = BlochWavefun(Qcell,n1,n2,n3,idxnz,mol.wks,mol.nspin);
        H.eband = cell(nkpts,1);
        for ik = 1:nkpts
            [X0{ik}, H.eband{ik}] = rotate_wfc(mol, H{ik}, X0{ik});
        end   
    else   
        X0 = Wavefun(Qcell{1},n1,n2,n3,idxnz);
        [X0, H.eband] = rotate_wfc(mol, H, X0); 
    end
end
end











