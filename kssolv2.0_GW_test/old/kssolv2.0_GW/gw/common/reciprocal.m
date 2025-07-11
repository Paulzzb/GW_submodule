function B2 = reciprocal(A)
% 输入：正晶格基矢量矩阵 A
% 输出：包含倒晶格基矢量的矩阵 B, 长度 lengths
%       B 矩阵的度规张量 bdot = B * B^T

%% 功能
% 1. 从正晶格到倒易晶格
% 2. 输入基矢量矩阵，并检查基矢量独立性

%% 未完成功能
% 1. 输入格式多样
% 2. 检查基矢量的线性独立性
noteString = '';
numNotPoints = 40;
for i = 1:numNotPoints
    noteString = [noteString '-'];
end

if (0)
% 从用户获取正晶格基矢量,并检查基矢量的独立性
    % 检查基矢量的线性独立性
    if det(A) < 3
        disp(exclamationString);
        disp('The base vector is dependent');
        disp(noteString)
    else
        disp('Direct vector')
        disp(A)
        disp(noteString)
        error('The base vector is independent'); % 如果矢量是线性独立的，退出循环
    end
end

% 1. 计算正晶格的体积
V = det(A);

% 2. 计算倒晶格基矢量矩阵，矢量定义和矩阵求逆,两种单位1/Angstron,2*pi/Angstron
%    A*B'=2*pi*I; B=2*pi*inv(A)'
%    注意 B 矩阵的转置
% % B = zeros(3, 3);
% % for i = 1:3
% %     B(i, :) = (2*pi/V) * cross(A(mod(i,3)+1, :), A(mod(i+1,3)+1, :));
% % end
B1 = 2*pi*inv(A)';
B2 = inv(A)'

% 使用norm函数计算每一行矢量的长度
lengths_B1 = zeros(1, size(B1, 1));
lengths_B2 = zeros(1, size(B2, 1));
for i = 1:size(B1, 1)
    lengths_B1(i) = norm(B1(i, :));
end
for i = 1:size(B2, 1)
    lengths_B2(i) = norm(B2(i, :));
end

% 3. 验证是否为对角矩阵 2*pi，这也是前面为什么可以求逆
disp(noteString)
disp("Please check it if it's a diagonal matrix")
disp("units:Angstron X 1/Angstron")
disp(B1*A')
disp(noteString)
disp("units:Angstron X 2*pi/Angstron")
disp(B2*A')
disp(noteString)

% 得到的矩阵 B 就是包含倒晶格基矢量的矩阵
disp(noteString)
disp('Reciprocal vectors (units:1/Angstron):')
disp(B1)
disp('The length of Rec-Vectors')
disp(lengths_B1)
disp(noteString)
disp("Reciprocal vectors (units:2*pi/Angstron):")
disp(B2)
disp('The length of Rec-Vectors')
disp(lengths_B2)
disp(noteString)
end