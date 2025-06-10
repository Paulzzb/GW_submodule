function config = read_input_param(filename)
% This function
% 1. Reads the input file from filename
% 2. Extracts the blocks and parameters
% 3. Check validation
% 4. Stores the parameters in a struct 'config'

% GET PARAM DEFINITION
def = allowed_param_list();

% READ FILE
fid = fopen(filename, 'r');
if fid == -1
  error('Cannot open file: %s', filename);
end
content = fread(fid, '*char')';
fclose(fid);

% Remove comments
content = regexprep(content, '%.*', '');

% Define pattern for blocks
blockPattern = '&(?<block>\w+)\s*\n(?<body>.*?)\nEND\s*&\k<block>';
% blockPattern = '&(?<block>\w+)[^&]*?END\s*&\k<block>';

% Extract blocks using case-insensitive match and dotall mode
matches = regexp(content, blockPattern, 'names');

% % Parse each match into block name and body
% blocks = struct([]);
% for i = 1:numel(matches)
%   % Extract block name
%   blkName = regexp(matches{i}, '&(\w+)', 'tokens', 'once');
%   body = regexp(matches{i}, '.*?\n(.*)\nEND\s*&\w+', 'tokens', 'once');
%   blocks(i).block = upper(blkName{1});
%   blocks(i).body = strtrim(body{1});
% end


% GET BLOCKS
% pattern = '&(?<block>\w+)(?<body>.*?)\/';
  

config = struct();
for i = 1:numel(matches)
  blockName = upper(matches(i).block);
  blockText = matches(i).body;

  if ~isfield(def, blockName)
  msg = sprintf('BLOCK NAME NOT DEFINED: %s', blockName);
  GWerror(msg);
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
    msg = sprintf('BLOCK "%s" HAS INVALID PARAMETER: %s', blockName, key);
    GWerror(msg);
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
