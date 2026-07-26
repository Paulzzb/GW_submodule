classdef SE_t < handle
    %SE_t  Mirror of Fortran QPN_m::SE_t (packages/QP/base/mod_QPN.F)
    %
    % Single-core / no-io use: set active=true, fill basis + PAR_* via
    % yambo_serial_par_scheme / helpers, then alloc_SF(NS_todo, NW).
    %
    % Not ported: full DESC_t, yMPI_comm, file I/O (io_SE).
    %
    % See also: yambo_serial_par_scheme, yambo_w_samp_static, example_init_cohsex_SE

    properties
        assigned (1,1) logical = false
        active   (1,1) logical = false
        kind     (1,:) char    = ''
        approx   (1,:) char    = ''
        basis    struct        % BASIS_t-like (see mod_BASIS.F)
        XC_id    (1,1) int32   = int32(0)
        NW       (1,1) int32   = int32(0)
        NG       (1,1) int32   = int32(0)
        NB       (1,2) int32   = int32([1, 1])   % band window for GW_driver
        desc     struct        % DESC_t placeholder (empty struct ok)
        prefactor (1,1) double = double(1)
        from_DB  (1,1) logical = false
        W        struct        % w_samp-like (mod_FREQUENCIES.F)
        SF                       % complex double, size [nState, nFreq] like Fortran SF(:,:)
        PAR_b    struct        % PAR_scheme serial stub
        PAR_qp   struct
        PAR_RL   struct
        PAR_q    struct
        PAR_structure (1,1) int32 = int32(0)
    end

    methods
        function obj = SE_t()
        end

        function alloc_SF(obj, nState, nFreq)
            validateattributes(nState, {'numeric'}, {'scalar','positive','integer'});
            validateattributes(nFreq, {'numeric'}, {'scalar','positive','integer'});
            obj.SF = complex(zeros(nState, nFreq, 'double'), zeros(nState, nFreq, 'double'));
        end

        function touch_assigned(obj)
            obj.assigned = true;
        end
    end
end
