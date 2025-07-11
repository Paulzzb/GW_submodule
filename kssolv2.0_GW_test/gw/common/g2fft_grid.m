function bidx = g2fft_grid(G_vector, n1, n2, n3)
    G_x = G_vector(:, 1);
    G_y = G_vector(:, 2);
    G_z = G_vector(:, 3);

    % 使用逻辑索引进行向量化操作
    G_x = modifyMatrixVectorized(G_x, n1);
    G_y = modifyMatrixVectorized(G_y, n2);
    G_z = modifyMatrixVectorized(G_z, n3);

    bidx = [G_x, G_y, G_z];
end

function result_matrix = modifyMatrixVectorized(matrix, x)
    % 使用逻辑索引进行向量化操作
    result_matrix = matrix + 1;  % 先对所有元素加1
    result_matrix(matrix < 0) = result_matrix(matrix < 0) + x;  % 对负数元素再加x
end