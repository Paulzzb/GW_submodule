function GWerror(msg)
% GWerror Display formatted error message and halt

    st = dbstack(1); % get caller info
    if isempty(st)
        loc = 'in anonymous context';
    else
        loc = sprintf('in %s at line %d', st(1).name, st(1).line);
    end

    fprintf(2, '[GW-ERROR] %s\n         --> %s\n', msg, loc);
    error('[GW] Execution terminated.');
end
