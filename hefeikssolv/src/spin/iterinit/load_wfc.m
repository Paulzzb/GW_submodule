function [mol, options]=load_wfc(type,qepath)
% only the initial wavefunction and ev are read
xmlname=[qepath,'/data-file-schema.xml'];
chargename=[qepath,'/charge-density.hdf5'];
data=xmlread(xmlname);
data=data.getElementsByTagName('output').item(0);

% ECUT
fft_grid=data.getElementsByTagName('fft_grid').item(0);
n1=str2double(fft_grid.getAttribute('nr1'));
n2=str2double(fft_grid.getAttribute('nr2'));
n3=str2double(fft_grid.getAttribute('nr3'));
ecutwfc=str2double(data.getElementsByTagName('ecutwfc').item(0).getTextContent);

% FUNCT
funct=data.getElementsByTagName('functional').item(0).getTextContent.toCharArray';
if strcmp(funct,'HSE')
    funct = 'HSE06'; % keep with kssolv
end

% nspin = 1/2/4, which represents calculation without spin, LSDA,
%  and noncolinear spin calculation.
lsda=data.getElementsByTagName('lsda').item(0).getTextContent;
noncolin=data.getElementsByTagName('noncolin').item(0).getTextContent;
lsda=strcmp(lsda,'true');
noncolin=strcmp(noncolin,'true');
% if soc is considerated 
spinorbit=data.getElementsByTagName('spinorbit').item(0).getTextContent;
spinorbit=strcmp(spinorbit,'true');
nspin = 1;
if lsda
    nspin = 2;
end
if noncolin
    nspin = 4;
end

% In the case of LSDA, nbnd_up = nbnd_dw
% Else only nbnd exists
if nspin == 2
    nbnd=str2double(data.getElementsByTagName('nbnd_up').item(0).getTextContent);
else
    nbnd=str2double(data.getElementsByTagName('nbnd').item(0).getTextContent);
end

% KPOINTS
% LSDA is a special case whose first half kpoints belong to up electrons
% and other belong to down electrons
nkpts=str2double(data.getElementsByTagName('nks').item(0).getTextContent);
kpoints_grid=data.getElementsByTagName('monkhorst_pack').item(0);
nk1=str2double(kpoints_grid.getAttribute('nk1'));
nk2=str2double(kpoints_grid.getAttribute('nk2'));
nk3=str2double(kpoints_grid.getAttribute('nk3'));

kptsobj=data.getElementsByTagName('k_point');
kpoints=zeros(nkpts,3);
for i=1:nkpts
    kpoints(i,:)=str2num(kptsobj.item(i-1).getTextContent); % unit is 2*pi/alat which is different from kssolv(a.u.)
end

% structure information
structure=data.getElementsByTagName('atomic_positions').item(0).getElementsByTagName('atom');
xyzlist=[];
atomlist=[];
for i=0:structure.getLength-1
    xyzlist=[xyzlist;str2num(structure.item(i).getTextContent)];
    atomlist=[atomlist, Atom(structure.item(i).getAttribute('name').toCharArray',0.5,90,0)]; % mag is not read from xml
end
% a1=Atom('Si',0);
% atomlist=[a1 a1 a1 a1 a1 a1 a1 a1];
cellobj=data.getElementsByTagName('cell').item(0);
a1=str2num(cellobj.getElementsByTagName('a1').item(0).getTextContent);
a2=str2num(cellobj.getElementsByTagName('a2').item(0).getTextContent);
a3=str2num(cellobj.getElementsByTagName('a3').item(0).getTextContent);
supercell=[a1;a2;a3];
atomstr=data.getElementsByTagName('atomic_structure').item(0);
alat=str2double(atomstr.getAttribute('alat')); % 2*pi/alat is the unit of reciprocal lattice vector

% electron temperature for broadening, fermi-dirac broadening is default
smear=data.getElementsByTagName('smearing').item(0);
if isempty(smear)
    degauss=0;
else
    degauss=str2double(smear.getAttribute('degauss'));
end
Temp = degauss*13.6/8.6173324*1e5;

% Build mol according to the value of type, Crystal or Molecular
if strcmp(type,'Molecule')
    mol = Molecule('supercell',supercell,'n1',n1,'n2',n2,'n3',n3,'atomlist',atomlist, 'xyzlist' ,xyzlist,...
        'ecut',ecutwfc,'funct',funct,'nbnd',nbnd,'nspin',nspin,'temperature',Temp,'lspinorb',spinorbit,'alat',alat);
elseif strcmp(type,'Crystal')
    mol = Crystal('supercell',supercell,'n1',n1,'n2',n2,'n3',n3,'atomlist',atomlist, 'xyzlist' ,xyzlist,...
        'ecut',ecutwfc,'funct',funct,'nbnd',nbnd,'nspin',nspin,'temperature',Temp,'lspinorb',spinorbit,...
        'autokpts',[nk1 nk2 nk3],'alat',alat);
end


if isempty(mol.ppvar)
    mol.ppvar  = PpVariable(mol);
    ne = 0;
     for it = 1:length(mol.atoms)
        ne = ne + mol.ppvar.venums(it)*mol.natoms(it);
    end
    if ne ~= mol.nel
        warning(['The number of valence electrons in Pp file is ' ...
            'different from Periodic Table']);
        mol = set(mol,'nel',ne);
    end
end

% RHO
% rho is divided into nspin parts
if nspin == 1
    millg=h5read(chargename,'/MillerIndices');
    millg=millg.';
    nl=mill2nl(millg,n1,n2,n3);
    rhog=h5read(chargename,'/rhotot_g');
    rhog=rhog(1:2:end-1)+1j*rhog(2:2:end);
    rhog3d=zeros(n1,n2,n3);
    rhog3d(nl)=rhog;
    rho=ifftn(rhog3d)*n1*n2*n3;
elseif nspin == 2
    rho = cell(2,1);
    % arho = rhoup + rhodw
    millg=h5read(chargename,'/MillerIndices');
    millg=millg.';
    nl=mill2nl(millg,n1,n2,n3);
    arhog=h5read(chargename,'/rhotot_g');
    arhog=arhog(1:2:end-1)+1j*arhog(2:2:end);
    arhog3d=zeros(n1,n2,n3);
    arhog3d(nl)=arhog;
    arho=real(ifftn(arhog3d)*n1*n2*n3);
    % drho = rhoup - rhodw
    drhog=h5read(chargename,'/rhodiff_g');
    drhog=drhog(1:2:end-1)+1j*drhog(2:2:end);
    drhog3d=zeros(n1,n2,n3);
    drhog3d(nl)=drhog;
    drho= real(ifftn(drhog3d)*n1*n2*n3);  
    %rho is stored by the format of up and down density
    rho{1} = (arho+drho)/2;
    rho{2} = (arho-drho)/2;
elseif nspin == 4
    rho = cell(4,1);
    %arho = rhoup + rhodw
    millg=h5read(chargename,'/MillerIndices');
    millg=millg.';
    nl=mill2nl(millg,n1,n2,n3);
    arhog=h5read(chargename,'/rhotot_g');
    arhog=arhog(1:2:end-1)+1j*arhog(2:2:end);
    arhog3d=zeros(n1,n2,n3);
    arhog3d(nl)=arhog;
    arho=real(ifftn(arhog3d)*n1*n2*n3);
    % mx
    %mxg=h5read(chargename,'/m_x');
    %mxg=mxg(1:2:end-1)+1j*mxg(2:2:end);
    %mxg3d=zeros(n1,n2,n3);
    %mxg3d(nl)=mxg;
    %mxr= real(ifftn(mxg3d)*n1*n2*n3);
    % my
    %myg=h5read(chargename,'/m_y');
    %myg=myg(1:2:end-1)+1j*myg(2:2:end);
    %myg3d=zeros(n1,n2,n3);
    %myg3d(nl)=myg;
    %myr= real(ifftn(myg3d)*n1*n2*n3);
    % mz
    %mzg=h5read(chargename,'/m_z');
    %mzg=mzg(1:2:end-1)+1j*mzg(2:2:end);
    %mzg3d=zeros(n1,n2,n3);
    %mzg3d(nl)=mzg;
    %mzr= real(ifftn(mzg3d)*n1*n2*n3);   
     rho{1} = arho;
    %rho{2} = mxr;
    %rho{3} = myr;
    %rho{4} = mzr;
     rho{2} = zeros(size(arho));
     rho{3} = zeros(size(arho));
     rho{4} = zeros(size(arho));
end

% EV DATA
ev = zeros(nbnd*(1+lsda),nkpts);
ev_data=data.getElementsByTagName('band_structure').item(0).getElementsByTagName('eigenvalues');
for ik = 1:nkpts
    ev(:,ik)=str2double(split(strtrim(string(ev_data.item(ik-1).getTextContent))));
end

% WAVEFUNCTION
% the truncation method is different in QE and kssolv, so we should 
% make a transformation(only reserve the values in idxnz of kssolv)
% when ecutwfc is large enough, the error will be small
grid = Ggrid(mol);
idxnz_ks = grid.idxnz;
npwx = length(idxnz_ks);
kpoints=kpoints*(2*pi/alat); % change unit to a.u.
if ~lsda  
    Qcell=cell(nkpts,1);
    for ik=1:nkpts
        wfcname=[qepath,'/wfc',num2str(ik),'.hdf5'];
        xk=h5readatt(wfcname,'/','xk')';
        assert(norm(xk-kpoints(ik,:))<1e-8,'Error: Kpoints in xml file and wavefunctions are different!')
        mill=h5read(wfcname,'/MillerIndices')';
        wfc=h5read(wfcname,'/evc');
        wfc=wfc(1:2:end-1,:)+1j*wfc(2:2:end,:);
        idxnz_qe=mill2nl(mill,n1,n2,n3);
        [~,iks,iqe]=intersect(idxnz_ks,idxnz_qe);
        wfctmp=zeros(npwx*(1+noncolin),nbnd);
        if ~noncolin
            wfctmp(iks,:)=wfc(iqe,:);
        else
            % noncolinear case is special for its wavefunctions can be
            % divided into up and down parts
            wfctmp(iks,:)=wfc(iqe,:);
            wfctmp(iks+npwx,:)=wfc(iqe+length(idxnz_qe),:);
        end   
        Qcell{ik}=wfctmp;
    end
else
    Qcell=cell(nkpts*2,1);
    % spin up
    for ik=1:nkpts
        wfcname=[qepath,'/wfcup',num2str(ik),'.hdf5'];
        xk=h5readatt(wfcname,'/','xk')';
        assert(norm(xk-kpoints(ik,:))<1e-8,'Error: Kpoints in xml file and wavefunctions are different!')
        mill=h5read(wfcname,'/MillerIndices')';
        wfc=h5read(wfcname,'/evc');
        wfc=wfc(1:2:end-1,:)+1j*wfc(2:2:end,:);
        idxnz_qe=mill2nl(mill,n1,n2,n3);
        [~,iks,iqe]=intersect(idxnz_ks,idxnz_qe);
        wfctmp=zeros(npwx,nbnd);
        wfctmp(iks,:)=wfc(iqe,:);
        Qcell{ik}=wfctmp;
    end
    % spin dw
    for ik=1:nkpts
        wfcname=[qepath,'/wfcdw',num2str(ik),'.hdf5'];
        xk=h5readatt(wfcname,'/','xk')';
        assert(norm(xk-kpoints(ik,:))<1e-8,'Error: Kpoints in xml file and wavefunctions are different!')
        mill=h5read(wfcname,'/MillerIndices')';
        wfc=h5read(wfcname,'/evc');
        wfc=wfc(1:2:end-1,:)+1j*wfc(2:2:end,:);
        idxnz_qe=mill2nl(mill,n1,n2,n3);
        [~,iks,iqe]=intersect(idxnz_ks,idxnz_qe);
        wfctmp=zeros(npwx,nbnd);
        wfctmp(iks,:)=wfc(iqe,:);
        Qcell{ik+nkpts}=wfctmp;
    end
end

% build wavefunctions and assign occupations
occ_data=data.getElementsByTagName('band_structure').item(0).getElementsByTagName('occupations');
if ~lsda
    if strcmp(type,'Molecule')
        X=Wavefun(Qcell{1},n1,n2,n3,idxnz_ks);
        X.occ=str2double(split(strtrim(string(occ_data.item(0).getTextContent))));
    elseif strcmp(type,'Crystal')
        X=BlochWavefun(Qcell,n1,n2,n3,idxnz_ks,mol.wks);
        for ik = 1:nkpts
            X{ik}.occ=str2double(split(strtrim(string(occ_data.item(ik-1).getTextContent))));
        end
    end
else
    if strcmp(type,'Molecule')
        X = cell(2,1);
        X{1} = Wavefun(Qcell{1},n1,n2,n3,idxnz_ks);
        X{2} = Wavefun(Qcell{2},n1,n2,n3,idxnz_ks);
        X{1} = set(X{1},'spin',1);
        X{2} = set(X{2},'spin',2);
        occ=str2double(split(strtrim(string(occ_data.item(0).getTextContent))));
        X{1}.occ=occ(1:nbnd);
        X{2}.occ=occ(nbnd+1:end);
    elseif strcmp(type,'Crystal')
        X = BlochWavefun(Qcell,n1,n2,n3,idxnz_ks,mol.wks,mol.nspin);
        for ik = 1:nkpts
            X{ik} = set(X{ik},'spin',1);
            X{ik+nkpts} = set(X{ik+nkpts},'spin',2);
            occ=str2double(split(strtrim(string(occ_data.item(ik-1).getTextContent))));
            X{ik}.occ=occ(1:nbnd);
            X{ik+nkpts}.occ=occ(nbnd+1:end);
        end        
    end
end

% set options
options = setksopt(mol);
% How to set rank for ISDF ???
% isdfoptions.rank=mol.nel*16*(1+noncolin)^2;
% isdfoptions.weight='power';
% isdfoptions.power=1;
% isdfoptions.init='random';
% isdfoptions.seed=0;
% isdfoptions.sys=mol;
% 
% options.isdfoptions=isdfoptions;
options.X0 = X;
options.ev0= ev;
options.rho0 = rho;
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
