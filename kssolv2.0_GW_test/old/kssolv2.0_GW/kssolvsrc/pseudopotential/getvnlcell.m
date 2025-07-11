function [vnlmatcell,vnlsigncell] = getvnlcell(cry)
% GETVNLCELL calculates the non-local pseudo potential for each k-points.
%    [vnlmatcell,vnlsigncell] = GETVNLCELL(cry,pseudovar) calculates the
%    non-local pseudo potential for each k-points in crystal.
%
%    See also VLOC2G, GETVNL.

nkpts = cry.nkpts;
vnlmatcell = cell(nkpts,1);
vnlsigncell = cell(nkpts,1);
for ik = 1:nkpts
    [vnlmatcell{ik},vnlsigncell{ik}] = ...
        getvnl(cry,cry.kpts(ik,:));
end

end
