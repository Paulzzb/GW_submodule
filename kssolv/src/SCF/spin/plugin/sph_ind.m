function sph_ind = sph_ind(l,j,m,spin)

if abs(j-l-0.5) < 1e-8
    if spin == 1
        sph_ind = m;
    elseif spin == 2
        sph_ind = m+1;
    end
elseif abs(j-l+0.5) < 1e-8
    if m < -l+1
        sph_ind = 0;
    else
        if spin == 1
            sph_ind = m-1;
        elseif spin == 2
            sph_ind = m;
        end
     end
end

if sph_ind<-l || sph_ind>l
    sph_ind = 0;
end
end

            
    