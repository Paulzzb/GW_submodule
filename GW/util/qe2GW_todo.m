function [GWinfor, sysinfo] = qe2GW(qepath)
% Set constants
ha2ry = 2.00;

sysinfo.vol       = [];
sysinfo.xyzlist   = [];
sysinfo.n1        = [];
sysinfo.n2        = [];
sysinfo.n3        = [];
sysinfo.ne        = [];
sysinfo.supercell = [];

GWinfor = GWinfo();

% filenames
xmlname=[qepath,'/data-file-schema.xml'];
chargename=[qepath,'/charge-density.hdf5'];
% Read vxc.dat
xml=[qepath,'/vxc.dat'];
fid = dlmread(xml);
[vxcrow,vxccol]=size(fid);
step=fid(1,4);
vxc.kpoints=[];
vxc.value =[];
for i =1:step+1:vxcrow-step
vxc.kpoints=[vxc.kpoints;fid(i,1),fid(i,2),fid(i,3)];
a=[];
  for j=i+1: 1 :i+step
      a=[a,fid(j,3)];
  end
vxc.value=[vxc.value;a];
end
vxc.value=vxc.value.';
GWinfor.Vxc = vxc;

% Read ev
ev=zeros(nbnd,nkpts);
ev_data=data.getElementsByTagName('band_structure').item(0).getElementsByTagName('eigenvalues');
for ik=1:nkpts
    ev(:,ik)=str2double(split(strtrim(string(ev_data.item(ik-1).getTextContent))));
end
GWinfor.ev = ev*ha2ry;


data=xmlread(xmlname);
data=data.getElementsByTagName('output').item(0);
fft_grid=data.getElementsByTagName('fft_grid').item(0);
% n1, n2, n3: the grid size in the real space.
n1 = str2double(fft_grid.getAttribute('nr1'));
n2 = str2double(fft_grid.getAttribute('nr2'));
n3 = str2double(fft_grid.getAttribute('nr3'));
% a1, a2, a3: the unit vectors in the real space.
cellobj=data.getElementsByTagName('cell').item(0);
a1 = str2num(cellobj.getElementsByTagName('a1').item(0).getTextContent);
a2 = str2num(cellobj.getElementsByTagName('a2').item(0).getTextContent);
a3 = str2num(cellobj.getElementsByTagName('a3').item(0).getTextContent);
supercell = [a1; a2; a3];
GWinfor.supercell = supercell;

% Ecut for scf-calculation
ecutwfc = str2double(data.getElementsByTagName('ecutwfc').item(0).getTextContent);
funct = data.getElementsByTagName('functional').item(0).getTextContent.toCharArray';
nbnd = str2double(data.getElementsByTagName('nbnd').item(0).getTextContent);
nqs = [];
nqs = [];
if strcmp(funct,'HSE')
  qgrid = data.getElementsByTagName('qpoint_grid').item(0).getAttributes;
  nqs = str2num([qgrid.item(0).getTextContent,qgrid.item(1).getTextContent,qgrid.item(2).getTextContent]);
end
nkpts = str2double(data.getElementsByTagName('nks').item(0).getTextContent);
%kptsobj = data.getElementsByTagName('starting_k_points').item(0).getElementsByTagName('k_point');
kptsobj = data.getElementsByTagName('k_point');
kpoints = zeros(nkpts, 3);
weight = zeros(nkpts, 1);
for i = 1:nkpts
  kpoints(i,:) = str2num(kptsobj.item(i-1).getTextContent);
  weight(i) = str2double(kptsobj.item(i-1).getAttributes.item(0).getTextContent);
end
weight = weight / sum(weight);
structure = data.getElementsByTagName('atomic_positions').item(0).getElementsByTagName('atom');
xyzlist = [];
atomlist = [];
for i = 0:structure.getLength-1
    xyzlist = [xyzlist; str2num(structure.item(i).getTextContent)];
    atomlist = [atomlist, Atom(structure.item(i).getAttribute('name').toCharArray')];
end


sys = Crystal('supercell',supercell,'n1',n1,'n2',n2,'n3',n3,'atomlist',atomlist, 'xyzlist' ,xyzlist,...
        'ecut',ecutwfc,'funct',funct,'nbnd',nbnd,'nkpts',nkpts,'kpts',kpoints,'wks',weight);

sysinfo.n1 = n1;
sysinfo.n2 = n2;
sysinfo.n3 = n3;
sysinfo.vol = det(supercell);
sysinfo.ne = 

% Read charge density information
ng = h5readatt(chargename,'/','ngm_g');
millg = h5read(chargename,'/MillerIndices');
millg = millg.';
nl = mill2nl(millg,n1,n2,n3);
rhog = h5read(chargename,'/rhotot_g');
rhog = rhog(1:2:end-1)+1j*rhog(2:2:end);
rhog3d = zeros(n1,n2,n3);
rhog3d(nl) = rhog;
rho = ifftn(rhog3d)*n1*n2*n3;
end
