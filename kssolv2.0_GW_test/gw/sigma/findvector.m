function [isortc] = findvector(G, gvec)
% This function calculates the index of a given reciprocal G in
% the space of E_cut.
nr = gvec.nr;
ng = gvec.ng;

index_vec = gvec.index_vec;
isortc = g2fft_index(G, gvec);
%isorti = zeros(length(G), 1);
del = zeros(length(G), 3);

% Vectorized conditions
valid_indices = (isortc >= 1 & isortc <= nr);
isortc(~valid_indices) = 0;
%isorti(~valid_indices) = 0;

% Update isortc with index_vec for valid indices
idx = isortc(valid_indices);
isortc(valid_indices) = index_vec(idx);

% Further filter valid indices within ng range
valid_ng_indices = (isortc >= 1 & isortc <= ng);
isortc(~valid_ng_indices) = 0;
%isorti(~valid_ng_indices) = 0;

% Update isorti for valid ng indices
%isorti(isortc(valid_ng_indices)) = find(valid_ng_indices);

% Calculate del for valid ng indices
G_valid = G(valid_ng_indices, :);
isortc_valid = isortc(valid_ng_indices);
mill_valid = gvec.mill(isortc_valid, :);
del(valid_ng_indices, :) = bsxfun(@minus, G_valid, mill_valid);

% Ensure G(i,:) == gvec.mill(isortc(i),:)
invalid_rows = any(del, 2);
isortc(invalid_rows) = 0;
%isorti(invalid_rows) = 0;
end