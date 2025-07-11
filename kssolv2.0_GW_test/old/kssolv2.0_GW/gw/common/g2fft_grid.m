function bidx = g2fft_grid(G_vector,n1,n2,n3)
G_x = G_vector(:,1);
G_y = G_vector(:,2);
G_z = G_vector(:,3);

G_x = modifyMatrix(G_x,n1);
G_y = modifyMatrix(G_y,n2);
G_z = modifyMatrix(G_z,n3);

bidx=[G_x G_y G_z];

function result_matrix = modifyMatrix(matrix,x)
    % 获取输入矩阵的大小
    [m, n] = size(matrix);
    
    % 初始化结果矩阵
    result_matrix = zeros(m, n);
    
    % 使用循环遍历输入矩阵的每个元素
    for i = 1:m
        for j = 1:n
            if matrix(i, j) < 0
                result_matrix(i, j) = matrix(i, j) + 1 + x;
            else
                result_matrix(i, j) = matrix(i, j) + 1;
            end
        end
    end
end
end