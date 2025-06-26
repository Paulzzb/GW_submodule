function display(obj)
% Custom display for GWOptions

    msg = ('--- GWOptions Summary ---\n');
    fprintf(msg);
    msg = sprintf('Still under development\n');
    QPerror(msg);
    % if isfield(obj.Constant, 'prefix')
        % fprintf('prefix: %s\n', obj.Constant.prefix);
    % end
    % if isfield(obj.GWCal, 'gw_method')
        % fprintf('GW method: %s\n', obj.GWCal.gw_method);
    % end
    % if isfield(obj.Groundstate, 'ecutwfc')
        % fprintf('ecutwfc: %.2f\n', obj.Groundstate.ecutwfc);
    % end
    % if isfield(obj.Groundstate, 'isGW')
        % fprintf('isGW: %d\n', obj.Groundstate.isGW);
    % end
end
