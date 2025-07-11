function [nrq,neq,indrk]=irrbz(syms,gr)
nrq=0;
fk=gr.f;
for ik=1:gr.nf
    found = 0;
    for it=1:syms.ntranq
        qk=gr.f(ik,:)*syms.mtrx{syms.indsub(it),:};
        [qk,kg]=krange(qk,1e-9);
        %% compare to other k-points in the irr. BZ with respect to qvec
        for irq=1:nrq
            if all(abs(fk(indrk(irq),:)-qk) < 1e-9)
                neq(irq) = neq(irq) + 1;
                found = 1;
                break
            end
        end
        if (found == 1)
            break
        end
    end
    if (found == 1)
        continue
    else
    nrq = nrq + 1;
    indrk(nrq) = ik; % vector connect idx of irr BZ and full BZ
    rq(nrq,:)=fk(ik,:);
    neq(nrq)=1;
    end
end
end