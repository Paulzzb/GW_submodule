function vcoul = generate_vcoul(Gvec, qvec, lattice, trunc_type, truncval)
% Generate v(q+G) Coulomb potential with truncation
% Inputs:
% - Gvec (ng x 3): integer reciprocal lattice vectors
% - qvec (1 x 3): q vector in Cartesian coordinates
% - lattice (3 x 3): real-space lattice vectors as rows
% - trunc_type: integer (0=no trunc, 4=wire, 5=box, 6=slab)
% - truncval: truncation parameters (e.g., length scale)
% Output:
% - vcoul: Coulomb potential values for each q+G

ng = size(Gvec,1);

% Compute reciprocal lattice
recip_lattice = 2*pi * inv(lattice');
Gcart = (recip_lattice * Gvec')';       % G vectors in Cartesian (ng x 3)
qG = Gcart + repmat(qvec, ng, 1);       % q+G vectors
qG2 = sum(qG.^2, 2);                    % |q+G|^2

% Initialize vcoul
vcoul = zeros(ng,1);

% Constants
fourpi = 4*pi;
eightpi = 8*pi;
cell_volume = abs(det(lattice));       % |a1 · (a2 x a3)|

% Default: unscreened Coulomb
switch trunc_type
    case 0  % No truncation (3D)
        vcoul = eightpi ./ qG2;        % v(q+G) = 8π / |q+G|²
        vcoul(qG2 == 0) = 0;           % avoid div-by-zero

    case 4  % Wire truncation (1D)
        % Apply wire-truncated Coulomb potential
        % Formula (e.g., from Rozzi et al. 2006):
        % v(q+G) = 4π / |q+G|² * [1 - exp(-|q+G|*L) * (|q+G|*L + 1)]
        L = truncval(1); % length along non-periodic directions
        qGnorm = sqrt(qG2);
        trunc_factor = 1 - exp(-qGnorm * L) .* (qGnorm * L + 1);
        vcoul = fourpi ./ qG2 .* trunc_factor;
        vcoul(qG2 == 0) = 0;

    case 5  % Box truncation (0D)
        % Use erfc-based cutoff (example form)
        L = truncval(1);
        qGnorm = sqrt(qG2);
        vcoul = fourpi ./ qG2 .* (1 - exp(-qGnorm.^2 * L^2));
        vcoul(qG2 == 0) = 0;

    case 6  % Slab truncation (2D)
        % Slab truncation along z (say Gz = G(:,3))
        Lz = truncval(3);
        Gz = qG(:,3);
        vcoul = eightpi ./ qG2 .* (1 - exp(-abs(Gz) * Lz));
        vcoul(qG2 == 0) = 0;

    otherwise
        error('Unsupported truncation type: %d', trunc_type);
end
end
