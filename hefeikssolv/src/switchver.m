function version = switchver(varargin)
%
% version = switchver('CPU'||'GPU')
% 
% Choose the kssolv version you want to run.The default version is CPU.
% We have transformed part of the Self Consistent Field (SCF) iteration code
% and the couplings code of TDDFT to GPU version.
% GPU codes include lobpcg, chebyfilt, davidson, davidson2, lanczos, omm

%
global version;
if nargin == 0
	version = 'CPU';
elseif varargin{1} == 'CPU'
	version = 'CPU';
elseif varargin{1} == 'GPU'
	version = 'GPU';
else
	error('Error. \nversion = switchver(''%s''). \nInput must be CPU or GPU or none!',varargin{1})
end
