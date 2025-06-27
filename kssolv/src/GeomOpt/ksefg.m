function [f, g, mol, info] = ksefg(x, varargin)

global saved_H saved_X relaxHistory;

nVarargs = length(varargin);
if nVarargs == 1
    mol = varargin{1}.mol;
    ksopts = varargin{1}.ksopts;
elseif nVarargs == 2
    mol = varargin{1};
    ksopts = varargin{2};
elseif nVarargs > 2
    error('Error: Too many input arguments')
end

natoms = sum(mol.natoms);
xyz = reshape(x, [natoms 3]);
mol = set(mol,'xyzlist',xyz);
if ~isempty(saved_H)
    ksopts.rho0 = saved_H.rho;
end
if ~isempty(saved_X)
    ksopts.X0 = saved_X;
end

[mol, saved_H, saved_X, info] = scf(mol, ksopts);
f = info.Etotvec(end);
xyzforces = mol.xyzforce;
g = -xyzforces(:);
fprintf('norm(g) = %11.3e\n', norm(g));

if size(relaxHistory,1) == 0
    relaxHistory(end+1,:) = {'Position','Energy','Force','Max F','RMS F','Displacement','Max D','RMS D',...
        'Matlab F','Matlab D','2-Norm F','2-Norm D','SCFIterNum'};
    relaxHistory(end+1,:) = {xyz,f,...
        ...% Tolerance used in Gaussian
        xyzforces,max(vecnorm(xyzforces,2,2)),norm(xyzforces(:))/sqrt(3*natoms),ones(natoms,3)*Inf,Inf,Inf,...
        ...% Tolerance used in MATLAB optimization toolbox
        norm(xyzforces(:),Inf),Inf,...
        ...% Tolerance used in other methods of KSSOLV
        norm(xyzforces(:)),Inf,...
        ...% SCF Iteration Number
        length(info.Etotvec)};
else
    deltaxyz = xyz - relaxHistory{end,1};
    relaxHistory(end+1,:) = {xyz,f,...
        ...% Tolerance used in Gaussian
        xyzforces,max(vecnorm(xyzforces,2,2)),norm(xyzforces(:))/sqrt(3*natoms),... % Maximum/RMS Force
        deltaxyz,max(vecnorm(deltaxyz,2,2)),norm(deltaxyz(:))/sqrt(3*natoms),...    % Maximum/RMS Displacement
        ...% Tolerance used in MATLAB optimization toolbox
        norm(xyzforces(:),Inf)/(1 + norm(relaxHistory{2,3}(:),Inf)),...             % OptimalityTolerance (TolFun)
        norm(deltaxyz(:)./(1+abs(xyz(:))),Inf),...                                  % StepTolerance (TolX)
        ...% Tolerance used in other methods of KSSOLV
        norm(xyzforces(:)),...                                                      % 2-Norm Force
        norm(deltaxyz(:)),...                                                       % 2-Norm Displacement
        ...% SCF Iteration Number
        length(info.Etotvec)};
end
%Threshold in Gaussian    default     tight
%Maximum Force           0.000450   0.000015
%RMS     Force           0.000300   0.000010
%Maximum Displacement    0.001800   0.000060
%RMS     Displacement    0.001200   0.000040
end