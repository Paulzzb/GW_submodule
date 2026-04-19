% License-Identifier: GPL
%
% Copyright (C) 2026
%
% Authors (see AUTHORS file for details): ZZ
%
% Last modified: 2026/03/18 ZZ

function getnGo()
  FFT_data = FFT.get();
  r_lat_data = lattice.manager('r_lat', 'get');
  
  indexGo = [];
  if ~isempty(r_lat_data.qindx_X)
    indexGo_tmp = ( r_lat_data.qindx_X(:, :, 2) );
    indexGo_tmp = indexGo_tmp(:);
    indexGo = unique([indexGo; indexGo_tmp]);
  end

  if ~isempty(r_lat_data.qindx_S)
    indexGo_tmp = ( r_lat_data.qindx_S(:, :, 2) );
    indexGo_tmp = indexGo_tmp(:);
    indexGo = unique([indexGo; indexGo_tmp]);
  end

  if ~isempty(r_lat_data.qindx_B)
    indexGo_tmp = ( r_lat_data.qindx_B(:, :, 2) );
    indexGo_tmp = indexGo_tmp(:);
    indexGo = unique([indexGo; indexGo_tmp]);
  end

  if ~isempty(r_lat_data.qindx_C)
    indexGo_tmp = ( r_lat_data.qindx_C(:, :, 2) );
    indexGo_tmp = indexGo_tmp(:);
    indexGo = unique([indexGo; indexGo_tmp]);
  end
  indexGo = sort(indexGo);

  FFT_data.nGo = max(indexGo);
  FFT_data.G_table = int32( zeros(r_lat_data.ng, FFT_data.nGo) );
  % FFT_data.indexGo = zeros(max(indexGo), 1, 'int32');
  %
  % for i = 1:FFT_data.nGo
  %   iGo = indexGo(i);
  %   FFT_data.indexGo(iGo) = int32(i);
  % end

  FFT.save2mod(FFT_data);

end
 