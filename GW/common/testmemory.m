
    [~, memInfo] = system('free -m');
    
    % 打印内存信息
    fprintf('当前内存信息:\n%s\n', memInfo);
    
    % 提取可用内存（可选）
    lines = strsplit(memInfo, newline);
    if length(lines) > 1
        availableMemoryLine = strsplit(lines{2});
        availableMemory = str2double(availableMemoryLine{7}); % 第七列是可用内存
        fprintf('当前可用内存: %.2f MB\n', availableMemory);
    end
