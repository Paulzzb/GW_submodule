cfile = mfilename('fullpath');
CPATH = fileparts(cfile);
CPATH = [CPATH, '/']; 
FILE_DIR = '../TMP_FILES/';

cd ../../
QPstartup
cd(CPATH);
%----- db test --------


root = fullfile(pwd, "ISDF_DB");
if ~exist(root, "dir"); mkdir(root); end


% ---- create meta skeleton (desc/xga/pmu/hVh + reserved G0set) ----
desc = struct();
desc.name = "ISDF database (skeleton)";
desc.created_at = char(datetime("now"));
desc.notes = "Fill in later";
meta = db_create(root, desc);


xga = rand(1000, 3, "double");          % example
pmu = complex(rand(2000, 64), rand(2000, 64)); % example complex
hVh = rand(64, 64, "single");           % example single

meta = db_write(root, meta, "xga", xga);
meta = db_write(root, meta, "pmu", pmu);
meta = db_write(root, meta, "hVh", hVh);

% ---- read back ----
xga2 = db_read(root, meta, "xga");
pmu2 = db_read(root, meta, "pmu");
hVh2 = db_read(root, meta, "hVh");

fprintf("xga diff: %.3e\n", norm(xga(:)-xga2(:)));
fprintf("pmu diff: %.3e\n", norm(pmu(:)-pmu2(:)));
fprintf("hVh diff: %.3e\n", norm(hVh(:)-hVh2(:)));

