function [Phi, Q] = scdm(Psi, varargin)

narginchk(1,3);
nargoutchk(1,2);

[N, nband] = size(Psi);

if nargin == 1
    [Q, R, piv] = qr(Psi',0);
    pivt(piv) = 1:N;
    Phi = R(:,pivt)';
end

% set parameters
switch nargin
    case 2
        tol = 5e-2;
    case 3
        tol = varargin{2};
end

if nargin > 1
    switch lower(varargin{1})
        case 'rand'
            rho = (Psi.^2)*ones(nband,1);
            [Q] = scdm_rand_transform(Psi,1e-8,nband,rho,1);
            Phi = Psi*Q;
        case 'fast'
            rho = (Psi.^2)*ones(nband,1);
            [Q_rand] = scdm_rand_transform(Psi,1e-8,nband,rho,1);
            Psi = Psi*Q_rand;
            [Q_fast] = scdm_fast_transform(Psi,1,tol);
            Phi = Psi*Q_fast;
            Q = Q_rand*Q_fast;
        otherwise
            error('Invalid algorithm choice, valid options are rand and fast\n');
    end
end
        