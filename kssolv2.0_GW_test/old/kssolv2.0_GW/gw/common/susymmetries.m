function umtrx = susymmetries(bvec, mtrx, itqq)
% This routine reads in a an SO(3) rotation matrix associated with a space group symmetry
% and generates the corresponding rotation matrix in the SU(2) group (in the postive branch).
% This matrix in SU(2), umtrx, is then used to rotate the spinors of the wavefunction.
% This routine is an essential part of all genwf routines, which use symmetry operations to
% reconstruct the wavefunctions in the little group of non-zero k- or q-vectors.
if itqq == 1
    umtrx = [1 0; 0 1];
else
    mtrx_tmp = mtrx;
    if det(mtrx_tmp) < 0 % for umtrx, improper rotations should be turned into proper rotations
        mtrx_tmp = -mtrx_tmp;
    end
    
    mtrx_cart = transpose(inv(bvec) * mtrx_tmp * bvec); % convert rotation matrix for lattice vectors to cartesian coordinates, R_cart =  B_Transpose * R * (B_Transpose)_Inverse
    
    [axis, angle] = rot_axis_angle(mtrx_cart);
    umtrx = su2_mat(axis, angle);
end

    function [axis, angle] = rot_axis_angle(mtrx)
        % Given a 3x3 orthogonal matrix, this routine gives the rotation axis and angle
        % This algorithm is based on the publication
        % F. Landis Markley, Journal of Guidance, Control, and Dynamics 31 (2008), p. 440
        % Ordering for quaternion: (axis(1) sin(angle/2), axis(2) sin(angle/2), axis(3) sin(angle/2), cos(angle/2))
        tr = trace(mtrx);
        cyclic = [2 3 1];
        xquat = zeros(4, 4);
        
        for ii = 1 : 4
            for jj = 1 : 4
                if (ii == jj && jj < 4)
                    xquat(ii, jj) = 1 + 2 * mtrx(ii, jj) - tr;
                elseif (ii ~= jj && ii < 4 && jj < 4)
                    xquat(ii, jj) = mtrx(ii, jj) + mtrx(jj, ii);
                elseif (ii == jj && jj == 4)
                    xquat(ii, jj) = 1 + tr;
                elseif (ii ~= jj && jj == 4)
                    xquat(ii, jj) = mtrx(6 - cyclic(ii) - ii, cyclic(ii)) - mtrx(cyclic(ii), 6 - cyclic(ii) - ii);
                end
            end
        end
        for ii = 1 : 4
            for jj = ii + 1 : 4
                xquat(jj, ii) = xquat(ii, jj);
            end
        end
        normx = sqrt(sum(xquat.^2, 2));
        maxnorm = max(normx); % Note that maxnorm is guaranteed to be at least 1.0d0
        indmax = normx == maxnorm;
        quat = xquat(:, indmax) / maxnorm;
        quat_pow = quat.^2;
        vecnorm = sqrt(sum(quat_pow(1:3)));
        angle = 2 * atan(vecnorm / quat(4));
        if vecnorm == 0
            axis = [0 0 1];
        else
            axis = quat(1 : 3) / vecnorm;
        end
    end

    function umtrx = su2_mat(axis, angle)
        % This routine takes in an axis and an angle from an SO(3) matrix and generates the corresponding SU(2) matrix.
        umtrx = zeros(2, 2);
        cosa = cos(angle * 0.5);
        sina = sin(angle * 0.5);
        umtrx(1, 1) = complex(cosa, -axis(3) * sina);
        umtrx(2, 1) = complex(-axis(2) * sina, -axis(1) * sina);
        umtrx(1, 2) = -conj(umtrx(2, 1));
        umtrx(2, 2) = conj(umtrx(1, 1)); % Note: this is actually the transpose of the usual U matrix, for use with genwf
    end
end