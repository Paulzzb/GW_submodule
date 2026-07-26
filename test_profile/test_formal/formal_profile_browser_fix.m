function index_path = formal_profile_browser_fix(profile_dir, matlab_root)
%FORMAL_PROFILE_BROWSER_FIX  Make profsave HTML viewable in a local browser.
%
%   formal_profile_browser_fix(profile_dir)
%
% Copies matlab-report-styles.css into profile_dir, rewrites file:// CSS links
% to a relative path, and writes index.html (alias of file0.html).

  if nargin < 1 || isempty(profile_dir)
    error('formal_profile_browser_fix:profile_dir', 'profile_dir is required.');
  end
  profile_dir = char(string(profile_dir));
  if ~isfolder(profile_dir)
    error('formal_profile_browser_fix:missing', 'Profile directory not found: %s', profile_dir);
  end

  if nargin < 2 || isempty(matlab_root)
    matlab_root = matlabroot;
  end
  css_name = 'matlab-report-styles.css';
  css_src = fullfile(matlab_root, 'toolbox', 'matlab', 'codetools', css_name);
  if ~isfile(css_src)
    error('formal_profile_browser_fix:css', 'MATLAB report CSS not found: %s', css_src);
  end
  css_dst = fullfile(profile_dir, css_name);
  if isfile(css_dst)
    delete(css_dst);
  end
  copyfile(css_src, css_dst, 'f');

  html_files = dir(fullfile(profile_dir, '*.html'));
  pat_abs = 'file:////public/software/MATLAB/R2023a/toolbox/matlab/codetools/matlab-report-styles.css';
  pat_rel = css_name;
  pat_generic = 'file://[^"]*matlab-report-styles\.css';

  for k = 1:numel(html_files)
    fpath = fullfile(profile_dir, html_files(k).name);
    txt = fileread(fpath);
    txt = regexprep(txt, pat_abs, pat_rel);
    txt = regexprep(txt, pat_generic, pat_rel);
    fid = fopen(fpath, 'w');
    if fid < 0
      error('formal_profile_browser_fix:write', 'Cannot write %s', fpath);
    end
    c = onCleanup(@() fclose(fid));
    fwrite(fid, txt);
    clear c
  end

  src0 = fullfile(profile_dir, 'file0.html');
  index_path = fullfile(profile_dir, 'index.html');
  if isfile(src0)
    copyfile(src0, index_path);
  end

  readme_path = fullfile(profile_dir, 'OPEN_IN_BROWSER.txt');
  fid = fopen(readme_path, 'w');
  if fid >= 0
    fprintf(fid, 'Open this folder in your local browser:\n\n');
    fprintf(fid, '  index.html          (profile summary)\n');
    fprintf(fid, '  file0.html          (same as index.html)\n');
    fprintf(fid, '  matlab-report-styles.css  (bundled stylesheet)\n\n');
    fprintf(fid, 'Remote server: download/copy this entire directory to your PC,\n');
    fprintf(fid, 'then double-click index.html, or run:\n');
    fprintf(fid, '  python3 -m http.server 8765\n');
    fprintf(fid, 'and visit http://localhost:8765/index.html\n');
    fclose(fid);
  end

  fprintf('formal_profile_browser_fix: %s -> open index.html\n', profile_dir);
end
