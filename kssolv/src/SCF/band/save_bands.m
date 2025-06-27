function save_bands(ebands,efermi,kpts,sys)
% save bands information to hdf5 file

fileID = H5F.create('KSSOLV_band.hdf5','H5F_ACC_TRUNC','H5P_DEFAULT','H5P_DEFAULT');

hdf5_double(fileID,ebands,'eigenvalues');
hdf5_attr_double(fileID,efermi,'efermi');
hdf5_double(fileID,sys.supercell,'cell');

symkpts=reshape([kpts{:,1:3}],[],3);
nkpts=[kpts{:,4}];
kpoints=[kpts{1,1:3}];
for i = 1:size(kpts,1)-1
    nk=nkpts(i)-1;
    dk=(symkpts(i+1,1:3)-symkpts(i,1:3)).*((1:nk)'/nk);
    kpoints=[kpoints;symkpts(i,1:3)+dk];
end
hdf5_double(fileID,kpoints,'kpoints');

labels=char(kpts{:,5});
hdf5_string(fileID,labels,'labels');
tick_idx=int32(cumsum([2,nkpts(1:end-1)]-1));
hdf5_int(fileID,tick_idx,'tick_idx');

H5F.close(fileID);
end

function hdf5_attr_double(fileID,attr,attr_name)
datatypeID= H5T.copy('H5T_NATIVE_DOUBLE');
dataspaceID = H5S.create('H5S_SCALAR');
attrID = H5A.create(fileID,attr_name,datatypeID,dataspaceID,'H5P_DEFAULT');
H5A.write(attrID,'H5ML_DEFAULT',attr)
end

function hdf5_double(fileID,array,array_name)
datatypeID = H5T.copy('H5T_NATIVE_DOUBLE');
dataspaceID = H5S.create_simple(ndims(array),size(array),[]);
datasetID = H5D.create(fileID,array_name,datatypeID,dataspaceID,'H5P_DEFAULT');
H5D.write(datasetID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT',array.');
end

function hdf5_int(fileID,array,array_name)
datatypeID = H5T.copy('H5T_NATIVE_INT32');
dataspaceID = H5S.create_simple(ndims(array),size(array),[]);
datasetID = H5D.create(fileID,array_name,datatypeID,dataspaceID,'H5P_DEFAULT');
H5D.write(datasetID,'H5ML_DEFAULT','H5S_ALL','H5S_ALL', 'H5P_DEFAULT',array.');
end

function hdf5_string(fileID,array,array_name)
datatypeID = H5T.copy ('H5T_C_S1');
H5T.set_size (datatypeID, size(array,2));
dataspaceID = H5S.create_simple(1, size(array,1), []);
datasetID = H5D.create (fileID, array_name, datatypeID, dataspaceID, 'H5P_DEFAULT');
H5D.write (datasetID, datatypeID, 'H5S_ALL', 'H5S_ALL', 'H5P_DEFAULT', array.');
end
