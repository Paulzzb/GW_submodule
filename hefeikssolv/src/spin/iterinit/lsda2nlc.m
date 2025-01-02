function flag = lsda2nlc(lsdapath,nlcpath)
% transform rho and X calculated from LSDA to the format suiting for
% noncolinear spin calculation

% RHO
chargename_lsda=[lsdapath,'/charge-density.hdf5'];
chargename_nlc=[nlcpath,'/charge-density.hdf5'];
arhog=h5read(chargename_lsda,'/rhotot_g');
drhog=h5read(chargename_lsda,'/rhodiff_g');
rho_xy=zeros(size(arhog));
h5write(chargename_nlc,'/rhotot_g',arhog);
h5write(chargename_nlc,'/m_x',rho_xy);
h5write(chargename_nlc,'/m_y',rho_xy);
h5write(chargename_nlc,'/m_z',drhog);
% WFC
wfcname=[lsdapath,'/wfcup1.hdf5'];
wfcup=h5read(wfcname,'/evc');
wfcname=[lsdapath,'/wfcdw1.hdf5'];
wfcdw=h5read(wfcname,'/evc');
wfcname=[nlcpath,'/wfc1.hdf5'];
npw = size(wfcup,1);
nbnd = size(wfcup,2);
count = size(wfcup);
wfc_updw = zeros(npw,nbnd);
h5write(wfcname,'/evc',wfcup,[1 1],count);
h5write(wfcname,'/evc',wfcdw,[npw+1 nbnd+1],count);
h5write(wfcname,'/evc',wfc_updw,[1 nbnd+1],count);
h5write(wfcname,'/evc',wfc_updw,[npw+1 1],count);
flag=1;
end