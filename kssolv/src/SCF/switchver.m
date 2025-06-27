function version = switchver(varargin)
% SWITCHVER Switch the version of KSSOLV to CPU or GPU
%    version = switchver('CPU'||'GPU') changes the global variable
%    version to CPU or GPU
% 
%    Choose the kssolv version you want to run.The default version is CPU.
%    We have transformed part of the Self Consistent Field (SCF) iteration code
%    and the couplings code of TDDFT to GPU version.
%    GPU codes include lobpcg, chebyfilt, davidson, davidson2, lanczos, omm

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

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
