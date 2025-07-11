cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/'];

cd ../../
QPstartup;
cd(CPATH)

input_driver('test');
fprintf('testqeread under construction\n');
