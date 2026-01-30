cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/']; 
FILE_DIR = '../TMP_FILES/';

cd ../../
QPstartup
cd(CPATH);
%----- db test --------

def = filename_map();

root = fullfile(pwd, def.isdf_database);
if ~exist(root, "dir"); mkdir(root); end


% ---- create meta skeleton (desc/xga/pmu/hVh + reserved G0set) ----
desc_d = desc();
desc_d.add("name", "ISDF database (skeleton)");
desc_d.add("notes", "Just a test");
% desc = struct();
% desc.name = "ISDF database (skeleton)";
% desc.created_at = char(datetime("now"));
% desc.notes = "Fill in later";
% desc.nv_range = 
meta = db_create(root, desc_d, 1);

xga = rand(1000, 3, "double");          % example
pmu = complex(rand(2000, 64), rand(2000, 64)); % example complex
hVh = rand(64, 64, "single");           % example single

ID = 1;
meta = db_write(root, meta, ID, "xga", xga);
meta = db_write(root, meta, ID, "pmu", pmu);
meta = db_write(root, meta, ID, "hVh", hVh);

% ---- read back ----
xga2 = db_read(root, ID, meta, "xga");
pmu2 = db_read(root, ID, meta, "pmu");
hVh2 = db_read(root, ID, meta, "hVh");

fprintf("xga diff: %.3e\n", norm(xga(:)-xga2(:)));
fprintf("pmu diff: %.3e\n", norm(pmu(:)-pmu2(:)));
fprintf("hVh diff: %.3e\n", norm(hVh(:)-hVh2(:)));

