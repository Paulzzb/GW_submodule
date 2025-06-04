function config = set_default_param_value(config, sys)
% This function fills in default values for missing config
% parameters with the help of system information
% Call after read_input_param

defaults = default_param_values();

blocks = fieldnames(defaults);
for i = 1:numel(blocks)
    blk = blocks{i};
    if ~isfield(config, blk)
        config.(blk) = struct();
    end
    keys = fieldnames(defaults.(blk));
    for j = 1:numel(keys)
        key = keys{j};
        if ~isfield(config.(blk), key)
            config.(blk).(key) = defaults.(blk).(key);
        end
    end
end