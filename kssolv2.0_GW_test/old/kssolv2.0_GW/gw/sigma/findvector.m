function[isortc,isorti]=findvector(G, gvec)
% This fuction can calculate the index of a given reciprical G in
% the space of E_cut.
n1=gvec.n1;
n2=gvec.n2;
n3=gvec.n3;
nr=gvec.nr;
ng=gvec.ng;

index_vec = gvec.index_vec;
isortc = g2fft_index(G, gvec);
isorti = zeros(length(G),1);
for i = 1:length(isortc)
    if isortc(i) >=1 && isortc(i) <= nr
        isortc(i) = index_vec(isortc(i));
        if isortc(i) >=1 && isortc(i) <= ng
            if (G(i,:) == gvec.mill(isortc(i),:))
                continue
            else
                isortc(i)=0;
            end
        else
            isortc(i)=0;
        end
    else
        isortc(i)=0;
    end
    isorti(isortc(i))=i;
end

