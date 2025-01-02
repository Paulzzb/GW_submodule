function ksinfo_new = setksinfo(ksinfo, options)
% set ksinfo with infos output by BerkeleyGW
% Those outputs are given by myself...to do the benchmark


% nameConstants = fieldnames(options.Constant);
% for i = 1:numel(nameConstants)
%     fieldname = nameConstants{i};
%     value = options.Constant.(fieldname);    
%     % Now, you can use the variable "fieldname" and "value" in your function
%     strEval = sprintf('%s = %f;', fieldname, value);
%     eval(strEval);
%     fprintf('Field: %s, Value: %s\n', fieldname, num2str(value));
% end
ry2ev = 13.605692530;
if ~isfield(options, 'setting_method')
	options.setting_method = 'kssolv';
end

nv = options.nv;
nc = options.nc;
ng_ksinfo = size(ksinfo.F);
for iv = 1 : nv+nc
  ksinfo.Z(:,iv) = ksinfo.Z(:,iv) / (norm(ksinfo.Z(:,iv)));
end

ksinfo_new = ksinfo;

switch lower(options.setting_method)
  case 'kssolv'
    if (1) % Copy wavefunctions from BGW for debugging purpose.
      gvec_bgw = readfile('gvec_components.csv', '%d %d %d');
      gvec_bgw = cell2mat(gvec_bgw);
      ks2bgw = index_ks2bgw(ksinfo.coulG(:, 1:3), gvec_bgw(1:ng_ksinfo, :));
      [~, bgw2ks] = sort(ks2bgw);
      
      wfnk.zk = readfile('wfnk.zk', '%f %f');
      wfnk.zk = reshape(complex(wfnk.zk{1}, wfnk.zk{2}), ng_ksinfo, []);
    	fprintf('Norm difference of Z = %f.\n', ...
    	        norm(wfnk.zk(bgw2ks, 1:options.nv+options.nc) - ksinfo.Z, 'fro'))
    	ksinfo.Z = wfnk.zk(bgw2ks, :);
      ev = readfile('wfnk.ek', '%f');
      ksinfo.ev = ev{1} / ry2ev;
    end
    
    gvec1 = ksinfo.gvec; 
    

    ksinfo_new.gvec = gvec1;
    ksinfo_new.Z = ksinfo.Z;
    ksinfo_new.ev = ksinfo.ev;

    % Z = ksinfo.Z * sqrt(ksinfo.vol);
    % for ind_aqs = 1:options.nv+options.nc
    % 	aqstmp = mtxel_sigma(ind_aqs, ksinfo, options);
    %   aqs_bgw = readfile(['aqs', num2str(ind_aqs)], '%f %f');
    % end
    % ksinfo_new.aqs = aqs;
    
	case 'bgw'
	  ksinfo_new = setksinfo_bgw(ksinfo, options);
  otherwise
		fprintf('options.setting_method = %s not support now!\n', options.setting_method)
	  error('Errors in setksinfo.')
	end


return
end

function ksinfo_new = setksinfo_bgw(ksinfo, options)
	ng_ksinfo = size(KSFFT(ksinfo.mol));


	ksinfo_new = ksinfo;
  rho = readfile('RHO.csv', '(%f, %f)');
  rho = rho{1} + sqrt(-1) * rho{2};
  vcoul = readfile('vcoul', '%f %f %f %d %d %d %f');
  % if (options.frequency_dependence == 1)
  gvec = struct();
  gvec.components = readmatrix('gvec_components.csv');
  index_vec  = readfile('gvec_index_vec.csv', '%d');
  gvec.index_vec  = index_vec{1};
  gvec.ng = length(find(gvec.index_vec));
  gvec.nfftgridpts = ksinfo.mol.n1 * ksinfo.mol.n2 * ksinfo.mol.n3;
  gvec.fftgrid = [ksinfo.mol.n1, ksinfo.mol.n2, ksinfo.mol.n3];
  % end 
  aqs_struct = dir([options.dir, '/aqs*']);
  aqs = cell(length(aqs_struct) - 1, 1);
  aqstmp = [];
  for ind_aqs = 1:8
    aqstmp = readfile(['aqs', num2str(ind_aqs)], '(%f, %f)');
    aqstmp = aqstmp{1} + sqrt(-1)*aqstmp{2};
    aqstmp = reshape(aqstmp, ng_ksinfo, []);
  	aqs{ind_aqs} = aqstmp;
  end
  
  ksinfo_new.rho  = rho;
	ksinfo_new.coulG = [double(vcoul{4}), double(vcoul{5}), double(vcoul{6}), vcoul{7}];  
  ksinfo_new.gvec = gvec;
  ksinfo_new.aqs = aqs;
end % function setksinfo_bgw
