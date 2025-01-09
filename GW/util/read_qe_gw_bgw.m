function [sys,options]=read_qe_gw_bgw(qepath)
% Output
%     sys: @Crystal
%     options: contains the following fields:
%              ev, rho0, X0(@Wavefunc), mill(G-grid)
xmlname=[qepath,'/data-file-schema.xml'];
chargename=[qepath,'/charge-density.hdf5'];
data=xmlread(xmlname);
data=data.getElementsByTagName('output').item(0);
fft_grid=data.getElementsByTagName('fft_grid').item(0);
n1=str2double(fft_grid.getAttribute('nr1'));
n2=str2double(fft_grid.getAttribute('nr2'));
n3=str2double(fft_grid.getAttribute('nr3'));
cellobj=data.getElementsByTagName('cell').item(0);
a1=str2num(cellobj.getElementsByTagName('a1').item(0).getTextContent);
a2=str2num(cellobj.getElementsByTagName('a2').item(0).getTextContent);
a3=str2num(cellobj.getElementsByTagName('a3').item(0).getTextContent);
supercell=[a1;a2;a3];
ecutwfc=str2double(data.getElementsByTagName('ecutwfc').item(0).getTextContent);
funct=data.getElementsByTagName('functional').item(0).getTextContent.toCharArray';
nbnd=str2double(data.getElementsByTagName('nbnd').item(0).getTextContent);
nqs=[];
nqs=[];
if strcmp(funct,'HSE')
    qgrid=data.getElementsByTagName('qpoint_grid').item(0).getAttributes;
    nqs=str2num([qgrid.item(0).getTextContent,qgrid.item(1).getTextContent,qgrid.item(2).getTextContent]);
end
nkpts=str2double(data.getElementsByTagName('nks').item(0).getTextContent);
%kptsobj=data.getElementsByTagName('starting_k_points').item(0).getElementsByTagName('k_point');
kptsobj=data.getElementsByTagName('k_point');
kpoints=zeros(nkpts,3);
weight=zeros(nkpts,1);
for i=1:nkpts
    kpoints(i,:)=str2num(kptsobj.item(i-1).getTextContent);
    weight(i)=str2double(kptsobj.item(i-1).getAttributes.item(0).getTextContent);
end
weight=weight/sum(weight);
structure=data.getElementsByTagName('atomic_positions').item(0).getElementsByTagName('atom');
xyzlist=[];
atomlist=[];
for i=0:structure.getLength-1
    xyzlist=[xyzlist;str2num(structure.item(i).getTextContent)];
    atomlist=[atomlist, Atom(structure.item(i).getAttribute('name').toCharArray')];
end
sys = Crystal('supercell',supercell,'n1',n1,'n2',n2,'n3',n3,'atomlist',atomlist, 'xyzlist' ,xyzlist,...
        'ecut',ecutwfc,'funct',funct,'nbnd',nbnd,'nkpts',nkpts,'kpts',kpoints,'wks',weight);
ev=zeros(nbnd,nkpts);
ev_data=data.getElementsByTagName('band_structure').item(0).getElementsByTagName('eigenvalues');
for ik=1:nkpts
    ev(:,ik)=str2double(split(strtrim(string(ev_data.item(ik-1).getTextContent))));
end
options.ev=ev*2;
ng=h5readatt(chargename,'/','ngm_g');
millg=h5read(chargename,'/MillerIndices');
millg=millg.';
nl=mill2nl(millg,n1,n2,n3);
rhog=h5read(chargename,'/rhotot_g');
rhog=rhog(1:2:end-1)+1j*rhog(2:2:end);
rhog3d=zeros(n1,n2,n3);
rhog3d(nl)=rhog;
rho=ifftn(rhog3d)*n1*n2*n3;



for ik=1:nkpts
    wfcname=[qepath,'/wfc',num2str(ik),'.hdf5'];
    millqe=h5read(wfcname,'/MillerIndices')';
    igwx=h5readatt(wfcname,'/','igwx');
    idxnzqe=zeros(igwx,1);
    idxnzqe(1:igwx, 1)=mill2nl(millqe,n1,n2,n3);%从小往大开始排，只取前igwx个有效地
    idxnz{1,ik}=idxnzqe;
    mill{1,ik}=millqe;
end
kpoints=kpoints/supercell'*2*pi;
for ik=1:nkpts
    wfcname=[qepath,'/wfc',num2str(ik),'.hdf5'];
    %h5disp(wfcname);
    xk=h5readatt(wfcname,'/','xk')';
    assert(norm(xk-kpoints(ik,:))<1e-8,'Error: Kpoints in xml file and wavefunctions are different!')
    mill=h5read(wfcname,'/MillerIndices')';
    nl=mill2nl(mill,n1,n2,n3);
    wfc=h5read(wfcname,'/evc');
    wfc=wfc(1:2:end-1,:)+1j*wfc(2:2:end,:);%波函数为复数
    %wfctmp=zeros(length(idxnz),size(wfc,2));
    %idx=map(nl);%按照n1里的顺序
    %wfctmp(idx,:)=wfc;%按照idxnz的顺序重新排
    %mill_k(idx,:)=mill;
    Qcell{ik}=wfc;
    mill_k{1,ik}=mill;
end

BX = BlochWavefun(Qcell,n1,n2,n3,idxnz,weight);
occ_data=data.getElementsByTagName('band_structure').item(0).getElementsByTagName('occupations');
for ik=1:nkpts
    BX{ik}.occ=str2double(split(strtrim(string(occ_data.item(ik-1).getTextContent))));
end

% HELLO!!!
options.idxnz = idxnz;
options.rho0 = rho;
options.X0 =BX;
options.mill=mill_k;

function nl=mill2nl(mill,n1,n2,n3)
%Convert mill index to nl (the index of the full G array)
assert(size(mill,2)==3,'Sencond dimension of mill should be 3!')
m1=mill(:,1);m2=mill(:,2);m3=mill(:,3);
m1=m1+int32((m1<0)*n1);
m2=m2+int32((m2<0)*n2);
m3=m3+int32((m3<0)*n3);
nl=m1+m2*n1+m3*n1*n2+1;
end

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
% Hello
% sys.vxc=vxc;

% GWoptions_in = [];
% GWinput = GWinfo();



end
