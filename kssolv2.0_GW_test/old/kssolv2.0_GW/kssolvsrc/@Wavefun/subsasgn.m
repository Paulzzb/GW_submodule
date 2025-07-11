function X = subsasgn(X,S,Y)
% WAVEFUN/SUBSASGN Subsasgn function for wave function class
%    X(idx) = Y assign the sub wave function.
%
%    See also Wavefun.

switch S(1).type
    case '()'
        rows = S(1).subs{1};
        cols = S(1).subs{2};
        if isa(Y,'Wavefun')
            X.psi(rows,cols) = Y.psi;
        else
            X.psi(rows,cols) = Y;
        end
    otherwise
        X = builtin('subsasgn',X,S,Y);
end

end
