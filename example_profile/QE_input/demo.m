cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/']; 

cd ../../
QPstartup
cd(CPATH);

input_driver('./test'); 
qp_driver('./SAVE'); % Input of qp_driver is the storage_dir in 'test'
