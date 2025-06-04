function config = read_input_param(filename)
% 读取输入文件并解析为结构体，同时进行参数合法性检查

% GET PARAM DEFINITION
def = param_def();

% READ FILE
fid = fopen(filename, 'r');
if fid == -1
    error('Cannot open file: %s', filename);
end
content = fread(fid, '*char')';
fclose(fid);

% GET BLOCKS
pattern = '&(?<block>\w+)(?<body>.*?)\/';
matches = regexp(content, pattern, 'names');

config = struct();
for i = 1:numel(matches)
    blockName = upper(matches(i).block);
    blockText = matches(i).body;

    if ~isfield(def, blockName)
        error('BLOCK NAME CANNOT BE IDENTIFIED: %s', blockName);
    end

    lines = regexp(blockText, '[\n\r]+', 'split');
    blockStruct = struct();

    for j = 1:numel(lines)
        line = strtrim(lines{j});
        if isempty(line) || startsWith(line, '!') || startsWith(line, '#')
            continue;
        end

        tokens = regexp(line, '(\w+)\s*=\s*(.*?)(,|$)', 'tokens');
        if isempty(tokens), continue; end

        key = lower(tokens{1}{1});
        valStr = strtrim(tokens{1}{2});

        if ~ismember(key, def.(blockName))
            error("BLOCK '%s' HAS INVALID PARAMETER: %s", blockName, key);
        end

        val = str2double(valStr);
        if isnan(val)
            if startsWith(valStr, '''') && endsWith(valStr, '''')
                val = valStr(2:end-1);
            else
                val = valStr;
            end
        end

        blockStruct.(key) = val;
    end

    config.(blockName) = blockStruct;
end
end
