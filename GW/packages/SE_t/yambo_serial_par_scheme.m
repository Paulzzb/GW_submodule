function scheme = yambo_serial_par_scheme(n_1d, n_table)
    %YAMBO_SERIAL_PAR_SCHEME  PAR_scheme stub for single CPU (mod_PARALLEL.F)
    %
    % Exposes fields read in packages/GW/GW_driver.F:
    %   scheme.IND.element_1D, scheme.IND.n_of_elements(scheme.ID+1)
    %   scheme.table(:)   for q mapping (identity 1:n_table if omitted)
    %   scheme.COMM_a2a.CPU_id == 0
    %
    % MPI communicators are not modeled; redux is no-op on one rank.

    if nargin < 2 || isempty(n_table)
        n_table = n_1d;
    end
    validateattributes(n_1d, {'numeric'}, {'scalar','positive','integer'});
    validateattributes(n_table, {'numeric'}, {'scalar','positive','integer'});

    scheme = struct();
    scheme.ID = int32(0);
    scheme.D = int32([1, 1]);
    scheme.N_ser = int32(n_1d);
    scheme.N_par = int32(1);
    scheme.COMM_world = int32(0);
    scheme.consecutive = true;
    scheme.IND = yambo_serial_pp_indexes(n_1d);
    scheme.table = int32(1:n_table);

    % yMPI_comm minimal stub (GW_driver checks COMM_a2a.CPU_id for I/O messaging)
    scheme.COMM_i = struct('CPU_id', int32(0), 'COMM', int32(0));
    scheme.COMM_a2a = struct('CPU_id', int32(0), 'COMM', int32(0));
end
