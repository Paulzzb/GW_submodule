% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/24

function out = manager(cmd, varargin)
	persistent wf_data

	if nargin == 0
		cmd = 'get';
	end

	switch lower(cmd)
		case 'get'
			if isempty(wf_data)
				error('wave_functions not initialized. Call wave_functions.save2mod or driver first.');
			end
			out = wf_data;
			return

		case 'save2mod'
			if isempty(varargin)
				error('wave_functions::save2mod requires one input argument.');
			end

			candidate = varargin{1};
			if ~(isa(candidate, 'wave_functions.base.wf_m'))
				error('wave_functions::save2mod requires input as a wave_functions.base.wf_m object.');
			end

			wf_data = candidate;
			out = int32(0);

		case 'free'
			wf_data = [];
			out = int32(0);

		otherwise
			error('Unknown command: %s', cmd);
	end
end
