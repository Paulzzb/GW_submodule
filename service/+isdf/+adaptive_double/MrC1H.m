function varargout = MrC1H(mode, varargin)
%MRC1H Persistent cache for MrC1H row blocks.
%   isdf.adaptive_double.MrC1H('ensure', Nw, Nisdfmax)
%   cols = isdf.adaptive_double.MrC1H('cached_cols', row_idx)
%   isdf.adaptive_double.MrC1H('set_range', row_idx, cstart, vals)
%   blk = isdf.adaptive_double.MrC1H('get', row_idx, cstart, cend)
%   isdf.adaptive_double.MrC1H('set_cached_cols', row_idx, ncols)
%   isdf.adaptive_double.MrC1H('clear')

  persistent cache cached_cols

  if nargin < 1
    error('MrC1H:mode', 'First argument ''mode'' is required.');
  end

  switch lower(mode)
    case 'ensure'
      if nargin < 3
        error('MrC1H:ensure', '''ensure'' requires Nw and Nisdfmax.');
      end
      Nw = double(varargin{1});
      Nisdfmax = double(varargin{2});
      if isempty(cache) || size(cache, 1) ~= Nw || size(cache, 2) ~= Nisdfmax
        cache = complex(zeros(Nw, Nisdfmax, 'double'));
        cached_cols = zeros(Nw, 1, 'uint32');
      end

    case 'cached_cols'
      row_idx = varargin{1}(:);
      varargout{1} = double(cached_cols(row_idx));

    case 'set_range'
      row_idx = varargin{1}(:);
      cstart = double(varargin{2});
      vals = varargin{3};
      cend = cstart + size(vals, 2) - 1;
      cache(row_idx, cstart:cend) = vals;

    case 'get'
      row_idx = varargin{1}(:);
      cstart = double(varargin{2});
      cend = double(varargin{3});
      varargout{1} = cache(row_idx, cstart:cend);

    case 'set_cached_cols'
      row_idx = varargin{1}(:);
      ncols = uint32(varargin{2});
      if isscalar(ncols)
        cached_cols(row_idx) = ncols;
      else
        cached_cols(row_idx) = ncols(:);
      end

    case 'clear'
      clear cache cached_cols

    otherwise
      error('MrC1H:mode', 'Unknown mode ''%s''.', mode);
  end
end
