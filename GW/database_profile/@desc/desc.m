% 
% License-Identifier: GPL
% 
% Copyright (C) 2026
% 
% Authors (see AUTHORS file for details): ZZ 
% 
classdef desc < handle
  %DESC Generic key-value descriptor with optional case-insensitive keys

  properties (Access = private)
    Map containers.Map
    CaseInsensitive logical = true 
  end

  methods
    function obj = desc(varargin)
      % Constructor
      %
      % Usage:
      %   d = desc();
      %   d = desc('CaseInsensitive', true);

      obj.Map = containers.Map('KeyType','char','ValueType','any');

      if nargin > 0
        if mod(nargin,2) ~= 0
          error('desc:Constructor', ...
            'Arguments must be name-value pairs.');
        end
        for i = 1:2:nargin
          name  = varargin{i};
          value = varargin{i+1};
          switch lower(string(name))
            case "caseinsensitive"
              obj.CaseInsensitive = logical(value);
            otherwise
              error('desc:Constructor', ...
                'Unknown parameter "%s".', name);
          end
        end
      end
    end

    function reset(obj)
      obj.Map = containers.Map('KeyType','char','ValueType','any');
    end

    function add(obj, name, value)
      key = obj.normalizeKey(name);
      obj.Map(key) = value;
    end

    function tf = has(obj, name)
      key = obj.normalizeKey(name);
      tf = isKey(obj.Map, key);
    end

    function value = get(obj, name, defaultValue)
      key = obj.normalizeKey(name);
      if isKey(obj.Map, key)
        value = obj.Map(key);
      else
        if nargin >= 3
          value = defaultValue;
        else
          error('desc:get:KeyNotFound', ...
            'Key "%s" not found.', key);
        end
      end
    end

    function remove(obj, name)
      key = obj.normalizeKey(name);
      if isKey(obj.Map, key)
        remove(obj.Map, key);
      end
    end

    function keysOut = keys(obj)
      keysOut = obj.Map.keys;
    end

    function n = count(obj)
      n = obj.Map.Count;
    end

    function tf = equals(obj, other)
      if ~isa(other, 'desc')
        tf = false; return;
      end

      % configuration must match
      if obj.CaseInsensitive ~= other.CaseInsensitive
        tf = false; return;
      end

      k1 = sort(obj.Map.keys);
      k2 = sort(other.Map.keys);

      if ~isequal(k1, k2)
        tf = false; return;
      end

      for i = 1:numel(k1)
        if ~desc.deepEqual(obj.Map(k1{i}), other.Map(k1{i}))
          tf = false; return;
        end
      end
      tf = true;
    end
  end

  methods (Access = private)
    function key = normalizeKey(obj, name)
      if isstring(name)
        if numel(name) ~= 1
          error('desc:KeyInvalid', ...
            'Key must be scalar string.');
        end
        key = char(name);
      elseif ischar(name)
        if size(name,1) ~= 1
          error('desc:KeyInvalid', ...
            'Key must be 1xN char.');
        end
        key = name;
      else
        error('desc:KeyInvalid', ...
          'Key must be string or char.');
      end

      if obj.CaseInsensitive
        key = lower(key);
      end
    end
  end

  methods (Static, Access = private)
    function tf = deepEqual(a, b)
      if isnumeric(a) || islogical(a)
        tf = isequaln(a, b); return;
      end

      if ischar(a) || isstring(a)
        tf = isequal(string(a), string(b)); return;
      end

      if iscell(a)
        if ~iscell(b) || ~isequal(size(a), size(b))
          tf = false; return;
        end
        tf = all(cellfun(@desc.deepEqual, a(:), b(:)));
        return;
      end

      if isstruct(a)
        if ~isstruct(b) || ~isequal(size(a), size(b))
          tf = false; return;
        end
        fa = sort(fieldnames(a));
        fb = sort(fieldnames(b));
        if ~isequal(fa, fb)
          tf = false; return;
        end
        for k = 1:numel(a)
          for i = 1:numel(fa)
            if ~desc.deepEqual(a(k).(fa{i}), b(k).(fa{i}))
              tf = false; return;
            end
          end
        end
        tf = true;
        return;
      end

      tf = isequaln(a, b);
    end
  end
end
