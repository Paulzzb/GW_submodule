
QPstartup

TMP_PATH = 'test_profile/TMP_FILES';
if ~exist(TMP_PATH, 'dir')
    mkdir(TMP_PATH);
end

cd test_profile/test_multik
test_input
% test_bz_samp
test_gw_x_k
cd ../../

% cd test_profile/test_groundstate 
% test_groundstate
% cd ../../

% cd test_profile/test_input
% test_input
% cd ../../

% cd test_profile/test_isdfdriver
% test_isdfdriver
% cd ../../

% cd test_profile/test_gw
% testgw
% cd ../../

% cd test_profile/testqpdir
% demo
% cd ../../

% cd test_profile/testqpISDF
% demo
% cd ../../

cd test