%kssolvpptype('pbe-hgh')
kssolvpptype('default')

%
%  construct bulk Silicon (crystal)
%
%  1. construct atoms
a1 = Atom('Si');
atomlist = [a1 a1 a1 a1 a1 a1 a1 a1];

%  2. set up the primitive cell
C = 10.216*eye(3);

%  3. define the coordinates the atoms
xyzlist = [
    0.1250    0.1250    0.1250
    0.1250    0.6250    0.6250
    0.6250    0.1250    0.6250
    0.6250    0.6250    0.1250
    0.8750    0.3750    0.8750
    0.3750    0.3750    0.3750
    0.3750    0.8750    0.8750
    0.8750    0.8750    0.3750
    ]*C;

% 4. Configure the molecule (crystal)
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',20,'name','Silicon');

opt = setksopt();
opt.maxscfiter = 30;

% QE ref: Etot = -63.48018027
[mol,H,X,info] = scf(mol,opt);
