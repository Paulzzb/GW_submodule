% One-shot: build SAVE/ from ./test (same pattern as test_Si/test_input.m).

cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
FILE_DIR = './SAVE/';

cd ../../
QPstartup
cd(CPATH);

service_reset_persistent();
packages_reset_persistent();
input_driver('./test');
load([FILE_DIR, 'GWinput.mat']);   % GWgroundstate
load([FILE_DIR, 'config.mat']);
