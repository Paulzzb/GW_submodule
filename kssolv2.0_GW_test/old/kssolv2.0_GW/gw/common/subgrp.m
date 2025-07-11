function [syms]=subgrp(qq, syms)
for it =1:syms.ntran
    qk = qq*syms.mtrx{it,1};
    dqk = qk - qq;
    [dqk, kg]=krange(dqk,1e-9);
    
    if (all(abs(dqk(:)) < 1e-9))
        syms.ntranq = syms.ntranq + 1;
        syms.indsub(syms.ntranq) = it;
        syms.kgzero(syms.ntranq,:) = -kg;
    end
end
end
