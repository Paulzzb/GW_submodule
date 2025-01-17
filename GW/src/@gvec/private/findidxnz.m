function idxnz = findidxnz(grid, n1n2n3)
% Find the mapping from 3-D grid indices to 1-D indices,
% Use 'our mapping rules'
n1 = n1n2n3(1);
n2 = n1n2n3(2);
n3 = n1n2n3(3);
grid = double(grid);


idxnz = 1 + (grid(:, 1) + (grid(:, 1) < 0) * n1) ...
        + (grid(:, 2) + (grid(:, 2) < 0) * n2) * n1 ... 
        + (grid(:, 3) + (grid(:, 3) < 0) * n3) * n1 * n2;

end % EOF