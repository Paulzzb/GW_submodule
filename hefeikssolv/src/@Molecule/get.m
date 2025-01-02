function val = get(a,attr_name)
%
% usage: val = get(a,attr_name);
%        retrieve various attributes of an Atom object a.
% e.g.  asymbol = get(a,'symbol') returns the atomic symbol associated with 
%       the Atom object a.
%       an = get(a,'anum') returns the atomic number associated with a
switch attr_name
  case {'name','Name'}
    val = a.name;
  case 'supercell'
    val = a.supercell;
  case 'alist'
    val = a.alist;
  case 'atoms'
    val = a.atoms;
  case 'xyzlist'
    val = a.xyzlist;
  case 'vol'
    val = a.vol;
  case 'natoms'
    val = a.natoms;
  case 'ecut'
    val = a.ecut;
  case 'ecut2'  
    val = a.ecut2;
  case 'n1'
    val = a.n1; 
  case 'n2'
    val = a.n2;
  case 'n3'
    val = a.n3;
  case 'rcut'
    val = a.rcut;
  case 'nel'
    val = a.nel;
  case 'vext'
    val = a.vext;
  case 'nspin'
    val = a.nspin;
  case 'k-points'
    val = a.kpts;
  case 'temperature'
    val = a.temperature;
  case 'funct'
    val = a.funct;
  case 'xyzforce'
    val = a.xyzforce;
  otherwise 
    error('invalid attribute requested');
end;

