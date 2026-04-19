function ik_map = map_kbz_by_symm(ikbz, isym, kptbz, rot_mtrx_RLU_G, tol)
S = double(rot_mtrx_RLU_G(:, :, isym));
target = double(kptbz(ikbz, :)) * S;

best = -1;
best_err = inf;
for ik = 1:size(kptbz, 1)
  diff = double(kptbz(ik, :)) - target;
  err = norm(diff - round(diff));
  if err < best_err
    best_err = err;
    best = ik;
  end
end

if best_err <= tol
  ik_map = int32(best);
else
  ik_map = int32(-1);
end
end
