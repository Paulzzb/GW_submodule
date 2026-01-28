% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZB 
%

function d = desc_add(d, name, value)
%DESC_add Add/overwrite (name,value) into a desc.
%
% Usage:
%   d = desc();
%   d = DESC_add(d, "version", 1.0);
%   d = DESC_add(d, "author", "Zhengbang Zhou");

    if nargin ~= 3
        error('DESC_add:Args', 'Usage: d = DESC_add(d, name, value)');
    end
    if ~isa(d, 'desc')
        error('DESC_add:Type', 'First input must be a desc.');
    end
    d.add(name, value);
end
