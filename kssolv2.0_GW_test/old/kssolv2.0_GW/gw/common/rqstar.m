function [nstar, indst, rqs]=rqstar(syms, rq)
nstar = 0;
for it = 1:syms.ntranq
    qk = rq * syms.mtrx{syms.indsub(it),:};
    [qk, gpt] = krange(qk, 1e-9);
    found = 0;
    for istar = 1:nstar
        if (all(abs(qk - rqs(istar, :)) < 1e-9)) 
            found = 1;
            break
        end
    end
    if (found == 0)
        nstar = nstar + 1;
        rqs(nstar, :) = qk;
        indst(nstar) = it;
    end
end
end