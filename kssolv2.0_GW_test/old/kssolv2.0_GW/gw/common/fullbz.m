function [gr]=fullbz(sys, syms, paranoid)
fk=zeros(syms.nrot,3);
gr.nr = sys.nkpts;
gr.r = round(sys.kpts / sys.bvec, 6); % reduce accuracy to prevent truncation errors
gr.nf=0;
for ir=1:gr.nr
    for it=1:syms.ntran
        tmpf=gr.r(ir,:)*syms.mtrx{it,1};
        [tmpf,gpt]=krange(tmpf,1e-9);
        found = 0;
        for ifull=1:gr.nf
            if all(abs(tmpf-fk(ifull,:)) < 1e-9)
                found = 1;
                break
            end
        end
    if (found == 1)
        continue
    else
        gr.nf=gr.nf+1;
        fk(gr.nf,:)=tmpf;
        itran(gr.nf)=it;
        indr(gr.nf)=ir;
    end
    end
end
gr.f=fk(1:gr.nf,:);
gr.itran=itran;
gr.indr=indr;
if (paranoid)
    for ii=1:gr.nf
        for jj=1:ii-1
            tmpf=abs(fk(ii,:)-fk(jj,:));
            tmpf=tmpf-floor(tmpf);
            tmpf(tmpf >= 0.5) = 1 - tmpf(tmpf >= 0.5);
            if sum(abs(tmpf)) <= 1e-9
                error('equivalent points found in the full BZ, equiv kpts %d and %d with diff %d', ii, jj, fk(ii,:)-fk(jj,:));
            end
        end
    end
end
end

