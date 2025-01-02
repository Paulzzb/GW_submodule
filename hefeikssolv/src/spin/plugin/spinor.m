function spinor = spinor(l,j,m,spin)

denom = 1/(2*l+1);
if abs(j-l-0.5) < 1e-8
    if spin == 1
        spinor = sqrt((l+m+1)*denom);
    elseif spin == 2
        spinor = sqrt((l-m)*denom);
    end
elseif abs(j-l+0.5) < 1e-8
    if m < -l+1
        spinor = 0;
    else
        if spin == 1
            spinor = sqrt((l-m+1)*denom);
        elseif spin == 2
            spinor = -sqrt((l+m)*denom);
        end
    end
end