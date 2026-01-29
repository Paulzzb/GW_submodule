
QPstartup

TMP_PATH = 'test_profile/TMP_FILES';
if ~exist(TMP_PATH, 'dir')
    mkdir(TMP_PATH);
end

cd test_profile/test_groundstate 
test_groundstate
cd ../../

cd test_profile/testInput
testinput
cd ../../

cd test_profile/test_isdfdriver
test_isdfdriver
cd ../../

cd test_profile/testgw
testgw
cd ../../
