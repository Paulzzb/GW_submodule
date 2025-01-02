function [Q, piv, conncomp, qr_size] = scdm_fast_transform(U,ortho,tol)

r = size(U,2);

Ua = U.^2;

col_cell = cell(r,1);
for k = 1:r
    Ua(:,k) = Ua(:,k) > tol*max(Ua(:,k));
    col_cell{k} = find(Ua(:,k));
end
Ua = sparse(double(Ua));
Dlogical = zeros(r,r);
Dlogical(Ua'*Ua > 0) = 1;
[perm, ~, conncomp] = dmperm(sparse(Dlogical));

piv_temp = [];
qr_size = zeros(length(conncomp)-1,2);

count_qr = 1;
for k = 1:(length(conncomp)-1)
    idx = conncomp(k):(conncomp(k+1)-1); 
    if length(idx) < 10
        cols = perm(idx);
        
        
        here = cell_union(col_cell(cols))';
        
        qr_size(count_qr,:) = [length(cols) length(here)];
        [~, ~, piv_loc] = qr(U(here,cols)',0);
        count_qr = count_qr+1;
        piv_temp = [piv_temp here(piv_loc(1:length(idx)))];
    else
        for j = idx
            cols = find(Dlogical(perm(j),:) > 0);
           
            
            here = cell_union(col_cell(cols))';
            
            qr_size(count_qr,:) = [length(cols), length(here)];
            [~, ~, piv_loc] = qr(U(here,cols)',0);
            count_qr = count_qr+1;
            piv_temp = [piv_temp here(piv_loc(1:length(cols)))];
        end
    end
end
[Q, ~, piv] = qr(U(piv_temp,:)',0);
piv = piv_temp(piv(1:r));
if ~ortho
    Q = U(piv,:)';
end