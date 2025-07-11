function [ikrkq, itqq, kgqq] = find_kpt_match(gr, syms, rkq)
found = 0;
for ik = 1:gr.nr
    for it = 1:syms.ntran
        qk = gr.r(ik, :) * syms.mtrx{it,:};
        %del = single(rkq - qk); %reduce accuracy
        del = rkq - qk;
        if (rem(del, 1) == 0) % determine whether the matrix elements are all integers
            found = 1;
            ikrkq = ik;
            itqq = it;
            kgqq = del;
            break
        end
        if (found == 1)
            break
        end
    end
end
        