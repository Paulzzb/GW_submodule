function [sys,options,syms,myneed]=read_qe_gw_bgw(qepath)
%The K point in xml and wfn differs by blat

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

lsda=data.getElementsByTagName('magnetization').item(0).getElementsByTagName('lsda').item(0).getTextContent;
noncolin=data.getElementsByTagName('magnetization').item(0).getElementsByTagName('noncolin').item(0).getTextContent;
if strcmpi(lsda, 'true')
    nspin = 2;
    if strcmpi(noncolin, 'true')
        nspinor = 2;
    else
        nspinor = 1;
    end
else
    nspin = 1;
    if strcmpi(noncolin, 'true')
        nspinor = 2;
    else
        nspinor = 1;
    end
end

nelec = str2double(data.getElementsByTagName('nelec').item(0).getTextContent);

if nspin == 2
    nbnd_up=str2double(data.getElementsByTagName('nbnd_up').item(0).getTextContent);
    nbnd_dw=str2double(data.getElementsByTagName('nbnd_dw').item(0).getTextContent);
    nbnd=(nbnd_up+nbnd_dw)/2;
else
    nbnd=str2double(data.getElementsByTagName('nbnd').item(0).getTextContent);
end

nqs=[];
if strcmp(funct,'HSE')
    qgrid=data.getElementsByTagName('qpoint_grid').item(0).getAttributes;
    nqs=str2num([qgrid.item(0).getTextContent,qgrid.item(1).getTextContent,qgrid.item(2).getTextContent]);
end
nkpts=str2double(data.getElementsByTagName('nk').item(0).getTextContent);
%kptsobj=data.getElementsByTagName('starting_k_points').item(0).getElementsByTagName('k_point');
kptsobj=data.getElementsByTagName('k_point');

%kpoints=zeros(nkpts,3);
weight=zeros(nkpts,1);
for i=1:nkpts
    %kpoints(i,:)=str2num(kptsobj.item(i-1).getTextContent);
    weight(i)=str2double(kptsobj.item(i-1).getAttributes.item(0).getTextContent);
end
weight=weight/sum(weight);
%kpoints = round(kpoints, 4);

structure=data.getElementsByTagName('atomic_positions').item(0).getElementsByTagName('atom');
xyzlist=[];
atomlist=[];
for i=0:structure.getLength-1
    xyzlist=[xyzlist;str2num(structure.item(i).getTextContent)];
    atomlist=[atomlist, Atom(structure.item(i).getAttribute('name').toCharArray')];
end

ev = zeros(nbnd,nkpts,nspin);

%%
nsym=str2double(data.getElementsByTagName('nsym').item(0).getTextContent);
nrot=str2double(data.getElementsByTagName('nrot').item(0).getTextContent);
rotation=data.getElementsByTagName('rotation');
mtrx=cell([nrot 1]);
for i=1:nrot
    mtrx{i,1}=str2num(rotation.item(i-1).getTextContent);
end

syms.ntran = nsym;
syms.ntranq = 0;
syms.mtrx = mtrx;
syms.nrot = nrot;
syms.indsub = zeros(nrot,1);
syms.kgzero = zeros(nrot,3);

%%

ev_data=data.getElementsByTagName('band_structure').item(0).getElementsByTagName('eigenvalues');

if nspin==2
    for ik=1:nkpts
        ev(:,ik,:)=reshape(str2double(split(strtrim(string(ev_data.item(ik-1).getTextContent)))),nbnd,[],2);
    end
else
    for ik=1:nkpts
        ev(:,ik)=str2double(split(strtrim(string(ev_data.item(ik-1).getTextContent))));
    end
end
options.ev = ev*2; %Unit from Ha to Ry

efermi = str2double(data.getElementsByTagName('fermi_energy').item(0).getTextContent);
options.efermi = efermi * 2; %Unit from Ha to Ry

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
    if nspin==2
        wfcname=[qepath,'/wfcdw',num2str(ik),'.hdf5'];
    else
        wfcname=[qepath,'/wfc',num2str(ik),'.hdf5'];
    end
    millqe=h5read(wfcname,'/MillerIndices')';
    igwx=h5readatt(wfcname,'/','igwx');
    idxnzqe=zeros(igwx,1);
    idxnzqe(1:igwx, 1)=mill2nl(millqe,n1,n2,n3);%从小往大开始排，只取前igwx个有效地
    idxnz{1,ik}=idxnzqe;
    mill{1,ik}=millqe;
end

Qcell=cell([nkpts nspin]);
kpoints = [];
kpoints_up = [];
kpoints_dw = [];

for ik=1:nkpts
    if nspin==2
        for ispin=1:nspin
            if ispin==1
                wfcname=[qepath,'/wfcup',num2str(ik),'.hdf5'];
                %h5disp(wfcname);
                xk=h5readatt(wfcname,'/','xk')';
                kpoints_up = [kpoints_up; xk];
                mill=h5read(wfcname,'/MillerIndices')';
                nl=mill2nl(mill,n1,n2,n3);
                wfc=h5read(wfcname,'/evc');
                wfc=wfc(1:2:end-1,:)+1j*wfc(2:2:end,:);%波函数为复数
                %wfctmp=zeros(length(idxnz),size(wfc,2));
                %idx=map(nl);%按照n1里的顺序
                %wfctmp(idx,:)=wfc;%按照idxnz的顺序重新排
                %mill_k(idx,:)=mill;
                
                Qcell{ik,ispin}=wfc;
                mill_k{1,ik}=mill;
            else
                wfcname=[qepath,'/wfcdw',num2str(ik),'.hdf5'];
                %h5disp(wfcname);
                xk=h5readatt(wfcname,'/','xk')';
                kpoints_dw = [kpoints_dw; xk];
                mill=h5read(wfcname,'/MillerIndices')';
                nl=mill2nl(mill,n1,n2,n3);
                wfc=h5read(wfcname,'/evc');
                wfc=wfc(1:2:end-1,:)+1j*wfc(2:2:end,:);%波函数为复数
                %wfctmp=zeros(length(idxnz),size(wfc,2));
                %idx=map(nl);%按照n1里的顺序
                %wfctmp(idx,:)=wfc;%按照idxnz的顺序重新排
                %mill_k(idx,:)=mill;
                
                Qcell{ik,ispin}=wfc;
                mill_k{1,ik}=mill;
            end
        end
        assert(isequal(kpoints_up, kpoints_dw),'The kpoints of spin up and spin down in WFN file are different!')
        kpoints = kpoints_up;
    else
        wfcname=[qepath,'/wfc',num2str(ik),'.hdf5'];
        %h5disp(wfcname);
        xk=h5readatt(wfcname,'/','xk')';
        kpoints = [kpoints; xk];
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
end

sys = Crystal('supercell',supercell,'n1',n1,'n2',n2,'n3',n3,'atomlist',atomlist, 'xyzlist' ,xyzlist,...
    'ecutwfc',ecutwfc,'funct',funct,'nbnd',nbnd,'nkpts',nkpts, 'kpts', kpoints, 'wks',weight,'nspin',nspin, 'nspinor', nspinor, 'nel', nelec);

BX = BlochWavefun(Qcell,n1,n2,n3,idxnz,weight);

occ_data=data.getElementsByTagName('band_structure').item(0).getElementsByTagName('occupations');

for ik = 1:nkpts
    occ = str2double(split(strtrim(string(occ_data.item(ik-1).getTextContent))));
    occ = reshape(occ, [], nspin);
    for ispin = 1:nspin
        ifmax = find(occ(:, ispin) >= 0.5, 1, 'last' );
        BX.wavefuncell{ik, ispin}.occ = occ(:, ispin);
        options.ifmax(ik, ispin) = ifmax;
    end
end

sys.idxnz=idxnz;
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
sys.vxc=vxc;

myneed = struct();

% System information 
myneed.supercell = supercell;
myneed.vol = det(supercell);
myneed.ne = nelec;
myneed.nspin = nspin;
myneed.nspinor = nspinor;
myneed.nkibz = nkpts;
myneed.n1 = n1;
myneed.n2 = n2;
myneed.n3 = n3;

% Things to recover the wavefunctions on real space.
myneed.wavefuncg = Qcell;
myneed.mill = mill_k;
myneed.idxnz = nl;

% kpoints
myneed.kibz = kpoints;
gr = fullbz(sys, syms, true);
myneed.kbz = gr.f;
myneed.ind_ikbz2ikibz = gr.indr;
myneed.indrot_ikbz2ikibz = gr.itran;
myneed.G0list_ikbz2ikibz = gr.gptlist;


% Things contains symmetric matrix information
% kpoints information


myneed.rho0 = rho; 
myneed.ev = ev*2;
vxc_ = zeros(nbnd,nkpts,nspin);
for ikpts = 1:nkpts
  if (nspin == 1)
    vxc_(:, ikpts) = vxc.value(1:nspin:end, ikpts);
  else
    for ispin = 1:nspin
      vxc_(:, ikpts, ispin) = vxc.value(ispin:nspin:end, ikpts);
    end
  end
end
myneed.vxc = vxc_;

occ_ = zeros(nbnd,nkpts,nspin);
for ik = 1:nkpts
    occ = str2double(split(strtrim(string(occ_data.item(ik-1).getTextContent))));
    occ = reshape(occ, [], nspin);
    occ_(:, ik, :) = occ;
    % for ispin = 1:nspin
    %     ifmax = find(occ(:, ispin) >= 0.5, 1, 'last' );
    %     BX.wavefuncell{ik, ispin}.occ = occ(:, ispin);
    %     options.ifmax(ik, ispin) = ifmax;
    % end
end
myneed.occ = occ_;
end
