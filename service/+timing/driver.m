% License-Identifier: BSD-3-Clause
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/04/08
%
% timing.driver() — allocate global/internal clock lists (Yambo timing_allocate).
% Does not read GW config; optional: timing.driver('nclock_max', N).
% MATLAB R2008a: inputParser.addParamValue (not addParameter).

function driver(varargin)
  p = inputParser;
  p.FunctionName = mfilename;
  p.addParamValue('nclock_max', [], @(x) isempty(x) || (isnumeric(x) && isequal(size(x), [1, 1]) && isfinite(x) && x > 0 && floor(x) == x));
  p.parse(varargin{:});
  nclock_max = p.Results.nclock_max;

  if isempty(nclock_max)
    tm = timing.base.timing_m.create_default();
  else
    tm = timing.base.timing_m.create_default(int32(nclock_max));
  end

  timing.save2mod(tm);
end
