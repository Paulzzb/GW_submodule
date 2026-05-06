function IND = yambo_serial_pp_indexes(n)
    %YAMBO_SERIAL_PP_INDEXES  Single-rank PP_indexes (mod_PARALLEL.F / PP_indexes)
    %
    % Fortran layout: n_of_elements(:) allocated to nCPU; only entry myid+1 used.
    % Serial: one rank, ID=0, everyone owns all indices.

    validateattributes(n, {'numeric'}, {'scalar','positive','integer'});

    IND = struct();
    IND.element_1D = true(1, n);
    IND.element_2D = true(n, n); %#ok<*NASGU> % optional superset
    IND.n_of_elements = int32(n);  % scalar: total work units (GW_driver uses one slot)
    IND.weight_1D = ones(1, n, 'int32');
    IND.first_of_1D = int32(1);
    IND.last_of_1D = int32(n);
end
