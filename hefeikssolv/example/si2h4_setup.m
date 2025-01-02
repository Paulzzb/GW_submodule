%
% construct the planar singlet silysilylene Si2H4 molecule (crystal)
%

%
% 1. construct atoms
%  J. Am. Chem. Soc.  p5222 (1979)

a1 = Atom('Si');
a2 = Atom('H');
atomlist = [a1 a1 a2 a2 a2 a2];

len_SiSi = 2.083;      % Si-Si bond length
len_SiH = 1.482;       % Si-H bond length
angle = 114.4/180*pi;  % H-Si-H angle

%
% 2. set up the supercell
%
BoxSize=10;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [      
-len_SiSi*0.5  0.0   0.0
 len_SiSi*0.5  0.0   0.0
 len_SiSi*0.5+len_SiH*cos(angle/2)  len_SiH*sin(angle/2)  0.0  
 len_SiSi*0.5+len_SiH*cos(angle/2) -len_SiH*sin(angle/2)  0.0
-len_SiSi*0.5-len_SiH*cos(angle/2)  len_SiH*sin(angle/2)  0.0  
-len_SiSi*0.5-len_SiH*cos(angle/2) -len_SiH*sin(angle/2)  0.0  
];
coefs = coefs/BoxSize;
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','Si2H4' );