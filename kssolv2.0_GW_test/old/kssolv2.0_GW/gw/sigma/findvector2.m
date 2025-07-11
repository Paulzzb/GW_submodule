function[isortc,isorti]=findvector2(G, gvec)
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
del = zeros(length(G),3);
for i = 1:length(isortc)
    if isortc(i) >=1 && isortc(i) <= nr
        isortc(i) = index_vec(isortc(i));
        if isortc(i) >=1 && isortc(i) <= ng
            isorti(isortc(i))=i;
            del(i, :) = G(i,:) - gvec.mill(isortc(i),:);
        else
            isortc(i)=0;
            isorti(i)=0;
        end
    else
        isortc(i)=0;
        isorti(i)=0;
    end
end
%% make sure G(i,:) == gvec.mill(isortc(i),:), if not:
[row, ~] = find(del);
isortc(row, :) = 0;
isorti(row, :) = 0;

