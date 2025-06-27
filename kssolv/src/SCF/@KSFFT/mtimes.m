function y = mtimes(F,x)
% MTIMES Overload multiplication operator for KSFFT class.
%    y = F*x returns the forward or inverse Fourier transform of x.
%
%    An important note !!!
%    ---------------------------------------------------------------------
%    In KSSOLV, the vector representing wave function in G-space
%    ( i.e. the coefficients of Fourier expansion ) is normalized:
%    	\sum_G |\psi(G)|^2 = 1
%
%    Therefore, to keep the integral normalization in real space, 
%    the wavefunction is writen as:
%	 	\psi(r) = 1/sqrt(V) \sum_G \psi(G) exp(iGr)
%    with
%    	\psi(G) = 1/sqrt(V) \int_\omega \psi(r) exp(-iGr) dr
%               = sqrt(V)/N \sum_r \psi(r) exp(-iGr)
%
%    The defination of FFT in MATLAB is:
%    	Y = fftn(X)  Y(k) = \sum_j X(j) W^(j-1)(k-1)_n
%       Y = ifftn(X) Y(j) = 1/N \sum_k Y(k) W^-(j-1)(k-1)_n
%       Wn = exp(-i*2pi/N)
%
%    Combining the above formulations, the transformation of wavefunction
%    between real space and G space can be calculated by:
%    	\psi(G) = sqrt(V)/N fftn[\psi(r)]
%       \psi(r) = N/sqrt(V) ifftn[\psi(G)]
%
%    The defination of KSFFT follows the above rules.
%    ---------------------------------------------------------------------
%
%    See also KSFFT.

%  Copyright (c) 2016-2017 Yingzhou Li and Chao Yang,
%                          Stanford University and Lawrence Berkeley
%                          National Laboratory
%  This file is distributed under the terms of the MIT License.

if (nargin == 2)
    idxnz = F.idxnz;
    n1 = F.n1;
    n2 = F.n2;
    n3 = F.n3;
    n123 = n1*n2*n3;
    vol = F.vol;
    if (isa(x,'numeric')||isa(x,'gpuArray'))
        [nrows, ncols ] = size(x);
        if (F.inverse)
            ng = length(idxnz);
            if ( nrows ~= ng )
                error('the number of rows in x does not match with the KSFFT object, nrows = %d, ng = %d', nrows, ng);
            end
            if (isa(x,'gpuArray'))
                a3d = gpuArray.zeros(n1,n2,n3);
                y = gpuArray.zeros(n1,n2,n3,ncols);
            else
                a3d = zeros(n1,n2,n3);
                y = zeros(n1,n2,n3,ncols);
            end
            for j = 1:ncols
                a3d(idxnz) = x(:,j);
                y(:,:,:,j) = ifftn(a3d);
                a3d(idxnz) = 0;
            end
            y = reshape(y,n123,ncols);
            y = y*n123/sqrt(vol);
        elseif (F.forward)
            if ( nrows ~= n123 )
                error('the number of rows in x does not match with the KSFFT object, nrows = %d', nrows, n123);
            end
            ng = length(idxnz);
            if (isa(x,'gpuArray'))
                y = gpuArray.zeros(ng,ncols);
            else
                y = zeros(ng,ncols);
            end
            for j = 1:ncols
                a3d = fftn(reshape(x(:,j),n1,n2,n3));
                y(:,j) = a3d(idxnz);
            end
            y = y*sqrt(vol)/n123;
        else
            error('KSFFT: something wrong with the FFT configuration');
        end
    else
        error('KSFFT must be applied to numeric data');
    end
else
    error('KSFFT syntax: y=F*x')
end
