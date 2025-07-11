function data = load_groundstate_info(dirin, typein)
% Load groundstate information and convert to @GWinfor format

% Step 1: Read raw data
switch lower(typein)
  case 'kssolv'
  data = load_kssolv_groundstate(dirin);
  case 'qe'
  data = load_qe_groundstate(dirin);
  otherwise
  msg = sprintf('Unsupported groundstate type: %s', typein);
  QPerror(msg);
end

end


function data = load_kssolv_groundstate(dirin)
% load_kssolv_groundstate - Load groundstate data from a KSSOLV-style .mat file
% User must have previously run a KSSOLV calculation and saved the results
% through save_groundstate_to_GWformat.m

filePath = fullfile(dirin, 'groundstate.mat');

if ~exist(filePath, 'file')
  msg = sprintf('groundstate.mat not found in directory: %s', dirin);
  QPerror(msg);
end

tmp = load(filePath);

if ~isfield(tmp, 'groundstate')
  msg = 'groundstate.mat does not contain a variable named "groundstate".';
  QPerror(msg);
end

required_fields = {'rhor', 'Vxc', 'ev', 'psig', 'sys', 'occupation', 'reciprocal_grid_info', 'nkpts', 'nspin', 'kpts', 'nspinor'};
for k = 1:length(required_fields)
  if ~isfield(tmp.groundstate, required_fields{k})
    msg = sprintf('Missing field "%s" in groundstate structure.', required_fields{k});
    QPerror(msg);
  end
end

data = tmp.groundstate;

end

function data = load_qe_groundstate(dirin)
% load_qe_groundstate - Load groundstate data from a QE pw.x output.
% User must have previously run a QE pw.x calculation.
myneed = load_qe_from_folder(dirin);
% error('load_qe_groundstate not implemented yet.');
data = struct();
data.rhor = myneed.rhor;
data.Vxc = myneed.vxc;
data.ev = myneed.ev;
data.psig = myneed.psig;
data.sys = myneed.sys;
fprintf("data.sys need construction.\n");
data.occupation = myneed.occupation;
reciprocal_grid_info = struct();
reciprocal_grid_info.fftgrid = [myneed.n1, myneed.n2, myneed.n3];
reciprocal_grid_info.vol = myneed.vol;
reciprocal_grid_info.idxnz = myneed.idxnz;
reciprocal_grid_info.wfncut = myneed.wfncut;
reciprocal_grid_info.xyz = myneed.mill;
fprintf("Later, check units of wfncut!\n");


data.reciprocal_grid_info = reciprocal_grid_info;


end

