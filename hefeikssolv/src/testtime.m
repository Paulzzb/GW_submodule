function time = testtime(varargin)
%
% version = switchver('CPU'||'GPU')
% 
% Choose the kssolv version you want to run.The default version is CPU.
% We have transformed part of the Self Consistent Field (SCF) iteration code
% and the couplings code of TDDFT to GPU version.
% GPU codes include lobpcg, chebyfilt, davidson, davidson2, lanczos, omm

%
global time;
if nargin == 0
	time = 'off';
elseif varargin{1} == 'on'
	time = 'on';
else
        time = 'off';
end
