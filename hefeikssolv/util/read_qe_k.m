function [sys,rho,BX,info]=read_qe_k(qepath,ng)
% Read in structure information, charge density and wave function of QE.
% This version is only for Crystal with k-points.
% In QE, idxnz varies with respect to k; while in KSSOLV, idxnz is invarient for different k.
% Here, idxnz is read and merged from QE, with some components padding with 0.
% Thus idxnz is different from the default value of KSSOLV, 'modifyidxnz'
% and 'qeidxnz' attributes are added, which should be used in FFT.
% ng=[ngx,ngy,ngz] is optional, it changes the default wavefunction grid size of QE.

xmlname=[qepath,'/data-file-schema.xml'];
chargename=[qepath,'/charge-density.hdf5'];
data=xmlread(xmlname).getElementsByTagName('output').item(0);
fft_grid=data.getElementsByTagName('fft_grid').item(0);
n1=str2double(fft_grid.getAttribute('nr1'));
n2=str2double(fft_grid.getAttribute('nr2'));
n3=str2double(fft_grid.getAttribute('nr3'));
if nargin>1
    n1=ng(1);
    n2=ng(2);
    n3=ng(3);
end
cellobj=data.getElementsByTagName('cell').item(0);
a1=str2num(cellobj.getElementsByTagName('a1').item(0).getTextContent);
a2=str2num(cellobj.getElementsByTagName('a2').item(0).getTextContent);
a3=str2num(cellobj.getElementsByTagName('a3').item(0).getTextContent);
supercell=[a1;a2;a3];
ecutwfc=str2double(data.getElementsByTagName('ecutwfc').item(0).getTextContent);
funct=data.getElementsByTagName('functional').item(0).getTextContent.toCharArray';
nbnd=str2double(data.getElementsByTagName('nbnd').item(0).getTextContent);
alat=str2double(data.getElementsByTagName('atomic_structure').item(0).getAttribute('alat'));
nqs=[];
if strcmp(funct,'HSE')
    qgrid=data.getElementsByTagName('qpoint_grid').item(0).getAttributes;
    nqs=str2num([qgrid.item(0).getTextContent,qgrid.item(1).getTextContent,qgrid.item(2).getTextContent]);
end
nkpts=str2double(data.getElementsByTagName('nks').item(0).getTextContent);
%kptsobj=data.getElementsByTagName('starting_k_points').item(0).getElementsByTagName('k_point');
startkptsobj=data.getElementsByTagName('monkhorst_pack').item(0);
nk1=str2double(startkptsobj.getAttribute('nk1'));
nk2=str2double(startkptsobj.getAttribute('nk2'));
nk3=str2double(startkptsobj.getAttribute('nk3'));
kptsobj=data.getElementsByTagName('k_point');
kpoints=zeros(nkpts,3);
weight=zeros(nkpts,1);
for i=1:nkpts
    kpoints(i,:)=str2num(kptsobj.item(i-1).getTextContent);
    weight(i)=str2double(kptsobj.item(i-1).getAttributes.item(0).getTextContent);
end
kpoints=kpoints/alat*supercell';
weight=weight/sum(weight);
structure=data.getElementsByTagName('atomic_positions').item(0).getElementsByTagName('atom');
xyzlist=[];
atomlist=[];
for i=0:structure.getLength-1
    xyzlist=[xyzlist;str2num(structure.item(i).getTextContent)];
    atomlist=[atomlist, Atom(structure.item(i).getAttribute('name').toCharArray')];
end

sys = Crystal('supercell',supercell,'n1',n1,'n2',n2,'n3',n3,'atomlist',atomlist, 'xyzlist' ,xyzlist,...
        'ecut',ecutwfc,'kpts',kpoints,'wks',weight,'funct',funct,'nqs',nqs,'nbnd',nbnd,'nkxyz',[nk1,nk2,nk3]);

ev=zeros(nbnd,nkpts);
ev_data=data.getElementsByTagName('band_structure').item(0).getElementsByTagName('eigenvalues');
for ik=1:nkpts
    ev(:,ik)=str2double(split(strtrim(string(ev_data.item(ik-1).getTextContent))));
end
info.ev=ev;
%info.force=str2num(data.getElementsByTagName('forces').item(0).getTextContent);

%h5disp(chargename);
ng=h5readatt(chargename,'/','ngm_g');
millg=h5read(chargename,'/MillerIndices');
millg=millg.';
nl=mill2nl(millg,n1,n2,n3);
rhog=h5read(chargename,'/rhotot_g');
rhog=rhog(1:2:end-1)+1j*rhog(2:2:end);
rhog3d=zeros(n1,n2,n3);
rhog3d(nl)=rhog;
rho=ifftn(rhog3d)*n1*n2*n3;

idxnz=zeros(ng,nkpts);
for ik=1:nkpts
    wfcname=[qepath,'/wfc',num2str(ik),'.hdf5'];
    mill=h5read(wfcname,'/MillerIndices')';
    igwx=h5readatt(wfcname,'/','igwx');
    idxnz(1:igwx, ik)=mill2nl(mill,n1,n2,n3);
end
idxnz=unique(idxnz);
idxnz=idxnz(2:end);
map=zeros(ng,1);
for i=1:length(idxnz)
    map(idxnz(i))=i;
end

Qcell = cell(nkpts,1);
for ik=1:nkpts
    wfcname=[qepath,'/wfc',num2str(ik),'.hdf5'];
    %h5disp(wfcname);
    xk=h5readatt(wfcname,'/','xk')';
    xk=xk*supercell'/2/pi;
    assert(norm(xk-kpoints(ik,:))<1e-8,'Error: Kpoints in xml file and wavefunctions are different!')
    mill=h5read(wfcname,'/MillerIndices')';
    nl=mill2nl(mill,n1,n2,n3);
    wfc=h5read(wfcname,'/evc');
    wfc=wfc(1:2:end-1,:)+1j*wfc(2:2:end,:);
    wfctmp=zeros(length(idxnz),size(wfc,2));
    idx=map(nl);
    wfctmp(idx,:)=wfc;
    Qcell{ik}=wfctmp;
end
BX = BlochWavefun(Qcell,n1,n2,n3,idxnz,weight);
occ_data=data.getElementsByTagName('band_structure').item(0).getElementsByTagName('occupations');
for ik=1:nkpts
    BX{ik}.occ=str2double(split(strtrim(string(occ_data.item(ik-1).getTextContent))));
end

sys.modifyidxnz=true;
sys.qeidxnz=idxnz;
end

function nl=mill2nl(mill,n1,n2,n3)
%Convert mill index to nl (the index of the full G array)
assert(size(mill,2)==3,'Sencond dimension of mill should be 3!')
m1=mill(:,1);m2=mill(:,2);m3=mill(:,3);
m1=m1+int32((m1<0)*n1);
m2=m2+int32((m2<0)*n2);
m3=m3+int32((m3<0)*n3);
nl=m1+m2*n1+m3*n1*n2+1;
end
