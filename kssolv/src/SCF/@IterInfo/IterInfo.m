classdef IterInfo
    % ITERINFO KSSOLV class for iteration information. 
    %    The IterInfo class contains the following fields.
    %        Field     Explaination
    %      ----------------------------------------------------------
    %        name        Crystal name
    %        converge    Indicator for convergence
    %        Eigvals     Eigenvalues at the final iteration
    %        Etotvec     Vector for total energy of each iteration
    %        SCFerrvec   Vector for SCF relative error of each iteration
    %        DCMerrvec   Vector for DCM relative error of each iteration
    %        TRDCMerrvec Vector for TRDCM relative error of each iteration
    %        Rhoerrvec   Vector for relative error of rho for each
    %                    iteration
    %        Vtoterrvec  Vector for relative error of Vtot for each
    %                    iteration
    %
    %    See also IterInfo.

    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    
    properties (SetAccess = public)
        converge = false
        Eigvals
        Etotvec
        SCFerrvec
        DCMerrvec
        TRDCMerrvec
    end
end
