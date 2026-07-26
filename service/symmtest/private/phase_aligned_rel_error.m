function rel_err = phase_aligned_rel_error(a, b)
a = double(a(:));
b = double(b(:));

if norm(a) == 0 && norm(b) == 0
  rel_err = 0.0;
  return;
end

alpha = sum(conj(a) .* b);
if abs(alpha) > 0
  b = b * exp(-1i * angle(alpha));
end

rel_err = norm(a - b) / max(norm(a), eps);
end
