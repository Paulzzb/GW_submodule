function y = mtimes(F,x)
%
% Perform the inverse Fourier transform from the reciprocal
% space to the real space
%
% usage: y = mtimes(F,x);
% The norm of operators:
% fft: sqrt(N)
% F: V/sqrt(N)
%
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
            y=reshape(y,n123,ncols);
            y = y*n123/vol;
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
            y = y*vol/n123;
        else
            error('KSFFT: something wrong with the FFT configuration');
        end
    else
        error('KSFFT must be applied to numeric data');
    end
else
    error('KSFFT syntax: y=F*x')
end
