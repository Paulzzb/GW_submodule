function config = read_input_param(filename)
    % 读取输入文件并解析为结构体
    fid = fopen(filename, 'r');
    if fid == -1
        error('无法打开文件: %s', filename);
    end
    content = fread(fid, '*char')';
    fclose(fid);

    % 查找所有 &BLOCK ... / 模块
    blockPattern = '&(?<block>\w+)(?<body>.*?)\/';
    matches = regexp(content, blockPattern, 'names');

    config = struct();
    for i = 1:numel(matches)
        blockName = upper(matches(i).block);
        blockText = matches(i).body;

        % 逐行解析 key = value
        lines = regexp(blockText, '[\n\r]+', 'split');
        blockStruct = struct();

        for j = 1:numel(lines)
            line = strtrim(lines{j});
            if isempty(line) || startsWith(line, '!') || startsWith(line, '#')
                continue;
            end

            % key = value[,]
            tokens = regexp(line, '(\w+)\s*=\s*(.*?)(,|$)', 'tokens');
            if isempty(tokens)
                continue;
            end

            key = tokens{1}{1};
            valStr = strtrim(tokens{1}{2});

            % 尝试转换为数字
            val = str2double(valStr);
            if isnan(val)
                % 去除引号（支持 'xxx'）
                if startsWith(valStr, '''') && endsWith(valStr, '''')
                    val = valStr(2:end-1);
                else
                    val = valStr; % 原样保留
                end
            end

            blockStruct.(key) = val;
        end

        config.(blockName) = blockStruct;
    end
end
