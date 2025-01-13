% test function for class gvec.
CPATH = mfilename('fullpath');
CPATH = fileparts(CPATH);
[GWPATH, ~, ~] = fileparts(CPATH);
CPATH = [CPATH, '/'];

cd(GWPATH); gw_startup(); cd(CPATH);

% Set Si8
n1 = 30; n2 = 30; n3 = 30;
supercell = [
   10.3340         0         0
         0   10.3340         0
         0         0   10.3340
];
ecut = 10.0;

gvecinput.ecut = ecut;
gvecinput.supercell = supercell;
gvecinput.n1 = n1;
gvecinput.n2 = n2;
gvecinput.n3 = n3;

gvec1 = gvec(gvecinput);

% Compare with a well-done one. 
GWinputfile = [CPATH, 'data/Si8/IntermediateFiles/GWinput_cohsex.mat'];
load(GWinputfile);
gvec1tocompare = GWinput.gvec;
propNames = properties(gvec1);

isEqual = true;
for i = 1:numel(propNames)
    if ~isequal(getfield(gvec1, propNames{i}), getfield(gvec1tocompare, propNames{i}))
        isEqual = false;
        disp(['gvec.' propNames{i} ' is not equal to the one in GWinput_cohsex.mat']);
        break;
    end
end


