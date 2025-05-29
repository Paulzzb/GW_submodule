function flag = issamegvec(gvec1, gvec2)
%   flag = issamegvec(gvec1, gvec2) returns true if gvec1 and gvec2 are the
%   same gvec.
%   The two gvecs are the same if they have the same
%   kpoint, same kpoint weight, and same basis set.

  % assert(isa(gvec1, 'gvec'));
  % assert(isa(gvec2, 'gvec'));

  TOL = 1e-8;
  flag = 1;


  fields1 = fieldnames(gvec1);
  fields2 = fieldnames(gvec2);
  if ~isequal(fields1, fields2)
    flag = 0;
    return;
  end 

  for i = 1:length(fields1)
    f = fields1{i};
    a = gvec1.(f);
    b = gvec2.(f);
    % Currently, all fields in gvec are numeric
    if isnumeric(a)
      if ~isequal(size(a), size(b))
        flag = 0;
        return;
      end
      if (norm(a(:) - b(:), 2) > TOL)
        flag = 0;
        return;
      end
    else
      error('Field %s is not comparable', f);
    end
  end

end