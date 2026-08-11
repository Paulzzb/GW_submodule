cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/']; 
FILE_DIR = './SAVE/';

cd ../../
QPstartup
cd(CPATH);

input_driver('./test');
load([FILE_DIR, 'GWinput.mat']);
GWinfo = GWgroundstate;
load([FILE_DIR, 'config.mat']);
