function degeneracy = degeneracy_detect()
  % Detect the degeneracy of a state specified by isc.
  % Return the number of degenerate states, including itself.

  system_data = system.get();
  nk = system_data.nk;
  nb = system_data.nb;

  %
  degen_tol = single(1e-5);
  

  first_idx_cell = cell(nk, 1);
  num_idx_cell = cell(nk, 1);
  degen_len = zeros(nk, 1, 'int32');
  for ik = 1:nk
    evals_k = system_data.Eo(:, ik, 1);
    [first_idx, num_idx] = build_degeneracy_segments(evals_k, degen_tol);
    first_idx_cell{ik} = first_idx;
    num_idx_cell{ik} = num_idx;
    degen_len(ik) = int32(numel(first_idx));
  end

  system_data.first_index_in_degeneracy = first_idx_cell;
  system_data.num_index_in_degeneracy = num_idx_cell;
  system_data.degeneracy_indices_len = degen_len;

  system.save2mod(system_data);
end

function [first_idx, num_idx] = build_degeneracy_segments(evals_k, tol)
  evals_k = single(evals_k(:));
  nb = numel(evals_k);

  first_idx = zeros(0, 1, 'int32');
  num_idx = zeros(0, 1, 'int32');
  if nb == 0
    return
  end

  seg_start = 1;
  for ib = 2:(nb + 1)
    segment_break = (ib == nb + 1);
    if ~segment_break
      segment_break = abs(evals_k(ib) - evals_k(ib - 1)) > tol;
    end

    if segment_break
      first_idx(end + 1, 1) = int32(seg_start); %#ok<AGROW>
      num_idx(end + 1, 1) = int32(ib - seg_start); %#ok<AGROW>
      seg_start = ib;
    end
  end
end