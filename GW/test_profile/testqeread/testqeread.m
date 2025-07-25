cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/'];

cd ../../
QPstartup;
cd(CPATH)

input_driver('test');
load ../TMP_FILES/GWinput.mat
fprintf('testqeread under construction\n');
