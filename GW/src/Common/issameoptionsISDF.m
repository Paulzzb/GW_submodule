function flag = issameoptionsISDF(type, op1, op2)
%   flag = issameoptionsISDF(op1, op2) returns true if op1 and op2 are the
%   same options.


  TOL = 1e-8;
  flag = 1;


  fields1 = fieldnames(op1);
  fields2 = fieldnames(op2);
  if ~isequal(fields1, fields2)
    flag = 0;
    return;
  end 

  for i = 1:length(fields1)
    f = fields1{i};
    a = op1.(f);
    b = op2.(f);
    % 跳过例外字段
    if strcmp(f, 'isdfoptions') || endsWith(f, 'Cauchy')
        continue;
    end
    switch type
      case 'vc'
        if (startsWith(f, 'vs') || startsWith(f, 'ss'))
          continue;
        end
      case 'vs'
        if (startsWith(f, 'vc') || startsWith(f, 'ss'))
          continue;
        end
      case 'ss'
        if (startsWith(f, 'vc') || startsWith(f, 'vs'))
          continue;
        end
    end
    % Currently, all fields in gvec are numeric
    if ischar(a) || isstring(a)
      if ~strcmp(a, b)
        flag = 0;
        return;
      end
    elseif isnumeric(a)
      if ~isequal(size(a), size(b))
        flag = 0;
        return;
      end
      if (norm(a(:) - b(:), 2) > TOL)
        flag = 0;
        return;
      end
    elseif islogical(a)
      if xor(a, b)
        flag = 0;
        return;
      end
    end
  end

end % function