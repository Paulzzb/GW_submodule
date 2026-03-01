cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/']; 
FILE_DIR = './TMP_FILES/';

cd ../../
QPstartup
cd(CPATH);

input_driver('./test');
load([FILE_DIR, 'GWinput.mat']);
isdf_driver(FILE_DIR)
qp_driver(FILE_DIR)
