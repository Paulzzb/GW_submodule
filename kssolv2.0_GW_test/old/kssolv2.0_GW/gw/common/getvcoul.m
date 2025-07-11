function coulG = getvcoul(ng, isrtq, ekin, amin, job_type)
gkk = ekin(isrtq(1 : ng));
for j = 1 : ng
    if ( abs(gkk(j)) ~= 0 )
        coulG(j, 1) = 8.0 * pi / (gkk(j));
        coulG(j, 1) = coulG(j, 1) * (1 - cos(sqrt(gkk(j)) * amin));
    else
        if job_type == 1
            coulG(j) = 8.0 * pi * amin^2 / 2;
        elseif job_type == 0
            coulG(j) = 0.0;
        end
    end
end

