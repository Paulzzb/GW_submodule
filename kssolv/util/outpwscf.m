function outpwscf(mol,filename)
%
% Purpose: Prepare an input file for PWSCF based on 
% the information contained in the Molecule object mol
%
% Usage: outpwscf(mol,filename);

molname = get(mol,'name');
alist   = get(mol,'alist');
xyzlist = get(mol,'xyzlist');
natoms  = size(xyzlist,1);
atoms  = get(mol,'atoms');
ntypes  = length(atoms);
C       = get(mol,'supercell');
fp = fopen(filename,'w');
for j = 1:110
   pseudoname{j} = 'unknown';
end;
pseudoname{1}  = 'H.vbc.UPF';
pseudoname{6}  = 'C.pz-vbc.UPF';
pseudoname{7}  = 'N.BLYP.UPF';
pseudoname{8}  = 'O.LDA.US.RRKJ3.UPF';
pseudoname{13} = 'Al.vbc.UPF';
pseudoname{14} = 'Si.vbc.UPF';
pseudoname{22} = 'Si.vdb.UPF';
pseudoname{33} = 'As.gon.UPF';
%
% write control statements and options
%
fprintf(fp,' &control\n');
fprintf(fp,'    calculation = ''scf''\n');
fprintf(fp,'    restart_mode= ''from_scratch'',\n');

fprintf(fp,'    prefix      = %s,\n', molname);
fprintf(fp,'    pseudo_dir  = ''/home/cyang/espresso-3.2.3/pseudo/''\n');
fprintf(fp,'    outdir      = ''/tmp''\n');
fprintf(fp,' /\n');
%
% write system statements and options
%
fprintf(fp,' &system\n');
fprintf(fp,'    ibrav = 0,\n');

fprintf(fp,'    nat   = %d,\n',natoms);
fprintf(fp,'    ntyp  = %d,\n',ntypes);
fprintf(fp,'    ecutwfc = 25.0,\n');
fprintf(fp,' /\n');
%
% write electrons statements and options
%
fprintf(fp,' &electrons\n');
fprintf(fp,'    diagonalization = ''cg''\n');
fprintf(fp,'    mixing_mode = ''plain''\n');
fprintf(fp,'    mixing_beta = 0.8\n');
fprintf(fp,'    conv_thr = 1.0d-8\n');
fprintf(fp,' /\n');
%
% write unit cell information
%
fprintf(fp,'CELL_PARAMETERS cubic\n');
for irow = 1:3
   fprintf(fp,' %12.4e   %12.4e   %12.4e\n', C(irow,1),C(irow,2),C(irow,3)); 
end;
%
% write atomic species and pseudo potential information
%
fprintf(fp,'ATOMIC_SPECIES\n');
for j = 1:ntypes
   anum = atoms(j).anum;
   symbol = atoms(j).symbol;
   amass = atoms(j).amass;
   fprintf(fp, ' %2s   %11.3e   %s\n', symbol, amass, pseudoname{anum}); 
end;
fprintf(fp,'ATOMIC_POSITIONS {bohr}\n');
for j = 1:natoms
   a = alist(j);
   xyz = xyzlist(j,:);
   symbol = atoms(a).symbol;
   fprintf(fp, ' %2s    %11.3e    %11.3e    %11.3e\n', ...
           symbol, xyz(1), xyz(2), xyz(3));
end;
fclose(fp);
