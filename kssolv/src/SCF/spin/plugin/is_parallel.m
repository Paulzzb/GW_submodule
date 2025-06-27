function is_parallel = is_parallel(u,v)
% judge if u and v are parallel.
cross = det([1 1 1;u';v']);
is_parallel = abs(cross)<1e-6;
end