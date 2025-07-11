function Y = subsref(X,S)
% WAVEFUN/SUBSREF Subsref function for wave function class
%    Y = X(idx) returns the sub wave function.
%
%    See also Wavefun.

if numel(S) > 1
    Y = builtin('subsref',X,S);
    return;
end

switch S(1).type
    case '()'
        rows = S(1).subs{1};
        cols = S(1).subs{2};
        Y = X;
        Y.psi = X.psi(rows,cols);
        if X.iscompact
            Y.idxnz = X.idxnz(rows);
        elseif length(rows) < Y.n1*Y.n2*Y.n3
            Y.iscompact = 1;
            Y.idxnz = rows;
        end
    otherwise
        Y = builtin('subsref',X,S);
end

end
