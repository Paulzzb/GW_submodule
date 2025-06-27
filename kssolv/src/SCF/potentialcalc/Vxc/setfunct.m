function factor = setfunct(funct)

switch upper(funct)
    case {'PZ','PW','PBE','SLA-PW-PBX-PBC'} % PZ: PRB 23, 5048 (1981); PW: PRB 45, 13244 (1992); PBE: PRL 77, 3865 (1996)
        factor.omega = 0;
        factor.exx_sr = 0;
        factor.exx_lr = 0;
        factor.dft_c = 1;
    case {'HSE06','HSE'} % J. Chem. Phys. 118, 8207–8215 (2003); J. Chem. Phys. 124, 219906 (2006); J. Chem. Phys. 125, 224106 (2006)
        factor.omega = 0.106; % 0.15/sqrt(2)[exx] and 0.15*2^(1/3)[dft_x] for HSE03; 0.11 for HSE06
        factor.exx_sr = 0.25;
        factor.exx_lr = 0;
        factor.dft_c = 1;
    case {'PBE0'} % J. Chem. Phys. 105, 9982–9985 (1996)
        factor.omega = 0;
        factor.exx_sr = 0.25;
        factor.exx_lr = 0.25;
        factor.dft_c = 1;
    case {'LC-WPBE'} % J. Chem. Phys. 125, 234109 (2006)
        factor.omega = 0.4;
        factor.exx_sr = 0;
        factor.exx_lr = 1;
        factor.dft_c = 1;
    case {'LC-PBE0'} % J. Chem. Phys. 129, 034107 (2008)
        factor.omega = 0.3;
        factor.exx_sr = 0.25;
        factor.exx_lr = 1;
        factor.dft_c = 1;
    case {'HF','PUREHF'}
        factor.omega = 0;
        factor.exx_sr = 1;
        factor.exx_lr = 1;
        factor.dft_c = 0;
    otherwise
        error('Unknown functional is used.')
end