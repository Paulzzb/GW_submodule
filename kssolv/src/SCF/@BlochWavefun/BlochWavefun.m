classdef BlochWavefun
    % BLOCHWAVEFUN KSSOLV class for Bloch wave function
    %    X = BLOCHWAVEFUN() returns an empty Bloch wave function.
    %
    %    X = BLOCHWAVEFUN(psi) constructs Bloch wave functions by cells of
    %    cells of 3D array psi. The first level of cell contains different
    %    wave functions whereas the second level contains the number of
    %    columns for each wave functions.
    %
    %    X = BLOCHWAVEFUN(psi,wks) constructs Bloch wave functions with wks
    %    being the weights for each k-points.
    %
    %    X = BLOCHWAVEFUN(psi,n1,n2,n3) returns nkpts Bloch wave functions
    %    and each of which contains ncols wave functions of size n1 by n2
    %    by n3 and filled by psi.
    %
    %    X = BLOCHWAVEFUN(psi,n1,n2,n3,wks) returns nkpts Bloch wave
    %    functions together with wks being the weights.
    %
    %    X = BLOCHWAVEFUN(psi,n1,n2,n3,idxnz) constructs compact Bloch wave
    %    functions by the psi with size n1 by n2 by n3. And idxnz indicates
    %    the indices for non-zero entries. Note the length of idxnz must
    %    not equal to the number of k-points.
    %
    %    X = BLOCHWAVEFUN(psi,n1,n2,n3,idxnz,wks) constructs compact Bloch
    %    wave functions with wks being the weights for k-points.
    %
    %    X = BLOCHWAVEFUN(psi,n1,n2,n3,idxnz,wks,nspin) constructs compact 
    %    Bloch wave functions with nspin being the number of density 
    %    components
    %
    %    Remark: This BlochWavefun serves as a data container for wave
    %    functions.
    %
    %    The Blochwavefun class contains the following fields.
    %        Field       Explaination
    %      ----------------------------------------------------------
    %        nspin       Number of density components
    %        nkpts       Number of k-points
    %        wks         The weight for each k-potins
    %        wavefuncell Cell of wave functions
    %
    %    See also Wavefun.
    
    %  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
    %                          Stanford University and Lawrence Berkeley
    %                          National Laboratory
    %  This file is distributed under the terms of the MIT License.
    properties (SetAccess = protected)
        nspin = 1
        nkpts = 0
        wks
        wavefuncell = {}
    end
    methods
        function BX = BlochWavefun(varargin)
            switch (nargin)
                case 0
                    return;
                case 1
                    Xin = varargin{1};
                    BX.nkpts = numel(Xin);
                    BX.wks = ones(BX.nkpts,1);
                    BX.wavefuncell = cellfun(@Wavefun,Xin, ...
                        'UniformOutput',0);
                case 2
                    Xin = varargin{1};
                    BX.nkpts = numel(Xin);
                    BX.wks = wks;
                    BX.wavefuncell = cellfun(@Wavefun,Xin, ...
                        'UniformOutput',0);
                case 4
                    Xin = varargin{1};
                    n1 = varargin{2};
                    n2 = varargin{3};
                    n3 = varargin{4};
                    BX.nkpts = numel(Xin);
                    BX.wks = ones(BX.nkpts,1);
                    BX.wavefuncell = ...
                        cellfun(@(X)Wavefun(X,n1,n2,n3),Xin, ...
                        'UniformOutput',0);
                case 5
                    Xin = varargin{1};
                    n1 = varargin{2};
                    n2 = varargin{3};
                    n3 = varargin{4};
                    vec = varargin{5};
                    BX.nkpts = numel(Xin);
                    if numel(vec) == numel(Xin)
                        BX.wks = vec;
                        BX.wavefuncell = ...
                            cellfun(@(X)Wavefun(X,n1,n2,n3),Xin, ...
                            'UniformOutput',0);
                    else
                        BX.wks = ones(BX.nkpts,1);
                        BX.wavefuncell = ...
                            cellfun(@(X)Wavefun(X,n1,n2,n3,vec),Xin, ...
                            'UniformOutput',0);
                    end
                case 6
                    Xin = varargin{1};
                    n1 = varargin{2};
                    n2 = varargin{3};
                    n3 = varargin{4};
                    idxnz = varargin{5};
                    BX.wks = varargin{6};
                    BX.nkpts = numel(Xin);
                    BX.wavefuncell = ...
                        cellfun(@(X)Wavefun(X,n1,n2,n3,idxnz),Xin, ...
                        'UniformOutput',0);
                case 7
                    Xin = varargin{1};
                    n1 = varargin{2};
                    n2 = varargin{3};
                    n3 = varargin{4};
                    idxnz = varargin{5};
                    BX.wks = varargin{6};
                    BX.nspin = varargin{7};
                    if BX.nspin ~= 2
                        BX.nkpts = numel(Xin);
                        BX.wavefuncell = ...
                        cellfun(@(X)Wavefun(X,n1,n2,n3,idxnz),Xin, ...
                        'UniformOutput',0);
                    else
                        nkpts = numel(Xin)/2;
                        BX.nkpts = nkpts;
                        X_up = cell(nkpts,1);
                        X_dw = cell(nkpts,1);
                        for ik = 1:nkpts
                            X_up{ik} = Xin{ik};
                            X_dw{ik} = Xin{ik+nkpts};
                        end
                        wavefuncell_up = ...
                        cellfun(@(X)Wavefun(X,n1,n2,n3,idxnz,1),X_up, ...
                        'UniformOutput',0);
                        wavefuncell_dw = ...
                        cellfun(@(X)Wavefun(X,n1,n2,n3,idxnz,2),X_dw, ...
                        'UniformOutput',0);
                        BX.wavefuncell = [wavefuncell_up;wavefuncell_dw];
                    end
                otherwise
                    error('Wrong number of arguements');
            end
        end
    end
end
