function flag = isdf_checkInputs(fName, type, Phir, nlist, mlist, gvec, vol, optionsISDF)
  TOL = 1e-8;
  flag = 1;

  % Construct struct
  Inputs = struct('Phir', Phir, 'nlist', nlist, 'mlist', mlist, 'gvec', gvec, 'vol', vol, 'optionsISDF', optionsISDF);
  % Check if file exists
  if isfile(fName)
    load(fName, 'ISDFinputs');
  else
    flag = 0;
    return
  end

  % If exists, check if inputs are the same
  fieldsA = fieldnames(Inputs);
  fieldsB = fieldnames(ISDFinputs);
  if ~isequal(fieldsA, fieldsB)
    flag = 0;
    return;
  end

  assert(isa(Inputs.gvec, 'gvec'));
  assert(isa(ISDFinputs.gvec, 'gvec'));
  % Compare the fields one by one
  for i = 1:length(fieldsA)
    f = fieldsA{i};
    a = Inputs.(f);
    b = ISDFinputs.(f);
    if isnumeric(a)
      if ~isequal(size(a), size(b))
        flag = 0;
        return;
      end
      if (norm(a(:) - b(:), 2) > TOL)
        flag = 0;
        return;
      end
    elseif isa(a, 'gvec')
      if ~issamegvec(a, b)
        flag = 0;
        return;
      end
    elseif strcmp(f, 'optionsISDF')
      if ~issameoptionsISDF(type, a, b)
        flag = 0;
        return;
      end
    else
      error('Field %s is not comparable', f);
    end
  end
end% function

