function basis = yambo_basis_minimal(NS_todo, table_nk, varargin)
    %YAMBO_BASIS_MINIMAL  BASIS_t-like struct (services/basis/base/mod_BASIS.F)
    %
    % NS_todo   number of QP transitions in this SE run
    % table_nk  int32 [NS_todo x 3] with columns [ib, ibp, ik_ibz] as in GW kinematics
    % Optional name-value:
    %   'NK', nk  (default size(table_nk,1) or 1)
    %
    % Fields used by GW_driver: NS_todo, states_todo, table

    p = inputParser;
    addRequired(p, 'NS_todo', @(x) isnumeric(x) && isscalar(x) && x > 0);
    addRequired(p, 'table_nk', @(x) isnumeric(x) && size(x, 2) == 3);
    addParameter(p, 'NK', [], @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    parse(p, NS_todo, table_nk, varargin{:});

    if int32(NS_todo) ~= int32(size(table_nk, 1))
        error('yambo_basis_minimal:NS_todo', 'NS_todo must equal size(table_nk,1).');
    end

    nk = p.Results.NK;
    if isempty(nk)
        nk = int32(max(table_nk(:, 3)));
        if nk < 1
            nk = int32(1);
        end
    else
        nk = int32(nk);
    end

    basis = struct();
    basis.reduction = double(100);
    basis.mixing = double(100);
    basis.E_range = double([-1, -1]);
    basis.mode = 'BZ';
    basis.IQ = int32(1);
    basis.NK = nk;
    basis.EH_pairs_only = true;
    basis.CB = int32([0, 0]);
    basis.VB = int32([0, 0]);
    basis.B = int32([0, 0]);
    basis.NB = int32(0);
    basis.NS = int32(NS_todo);
    basis.NS_serial = int32(NS_todo);
    basis.NS_todo = int32(NS_todo);
    basis.states_todo = int32(1:NS_todo);
    basis.table = int32(table_nk);
    basis.pos = [];
    basis.E = [];
end
