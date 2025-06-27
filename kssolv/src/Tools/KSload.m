function [mol,H,X,info] = KSload(filename)

ksdat = load(filename);
mol = ksdat.mol;
H = ksdat.H;
X = ksdat.X;
info = ksdat.info;

end