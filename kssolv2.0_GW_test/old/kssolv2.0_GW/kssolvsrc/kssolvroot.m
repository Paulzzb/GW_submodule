function path = kssolvroot()
% KSSOLVROOT Root directory of KSSOLV.
%   S = KSSOLVROOT returns a string that is the absolute path of the
%   directory where the KSSOLVE software is installed.
%  
%   See also matlabroot.

path = mfilename('fullpath');
path = path(1:end-20);

end
