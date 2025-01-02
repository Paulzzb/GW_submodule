clc;
clear all;
close all;

% Set timer
tstart  = cputime;
%
% construct the H2 molecule 
%
kssolvpptype('default');
%
% 1. construct atoms
%
a1 = Atom('H');
atomlist = [a1 a1];
%
% 2. set up supercell
%
C = 10*eye(3);
%
% 3. define the coordinates the atoms 
%
xyzlist = [
 1.5     0.0      0.0
 0.0     0.0      0.0
];

%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
               'ecut',40.0, 'name','h2' );

% 
% Set up MD run
% 
nAtoms = sum(mol.natoms);
ksopts = setksopt();
ksopts.maxscfiter = 200;
ksopts.betamix = 0.08;
ksopts.mixdim = 8;
ksopts.scftol = 1.0E-6;
ksopts = setksopt(ksopts,'mixtype','pulay','eigmethod','eigs');

% First SCF cycle:
[mol, H0, X0, info] = scf(mol, ksopts);

nMdSteps = 500;

[~, ~, ~, eTot, eNKin, INFO] = md(mol, ksopts, H0, X0, nMdSteps);

% [X, V, F, eTot, eKin, eFree, entropy, INFO] = md(mol, H0, X0, nMdSteps, ksopts);

% dlmcell('H2_MD_runSummary.txt', INFO);

% https://dgleich.wordpress.com/2013/06/04/creating-high-quality-graphics
% -in-matlab-for-papers-and-presentations/
% Defaults for this blog post
width = 6;     % Width in inches
height = 4;    % Height in inches
alw = 5.0;     % AxesLineWidth
fsz = 18;      % Fontsize
lw = 2.0;      % LineWidth
% 
% Plot Total energy:
% 
figure
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*1000, height*1000]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', fsz, 'FontWeight', 'Bold');
hold on;
% 
plot(1:nMdSteps, eTot/nAtoms, 'r', 'LineWidth', 1.25)
% 
xlabel('Time (fs)','Interpreter','tex');
ylabel('Total Energy (Ha/atom)','Interpreter','tex');
box on; 
axis square;
set(gca,'LineWidth',lw);
set(gca,'FontSize',fsz);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');

% 
% Plot Kinetic energy:
% 
figure
set(gcf,'InvertHardcopy','on');
set(gcf,'PaperUnits', 'inches');
pos = get(gcf, 'Position');
set(gcf, 'Position', [pos(1) pos(2) width*1000, height*1000]);
set(gca, 'FontSize', fsz, 'LineWidth', alw);
set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
set(get(gca,'title'),'FontSize', fsz, 'FontWeight', 'Bold');
hold on;
% 
plot(1:nMdSteps, eNKin/nAtoms, 'b', 'LineWidth', 1.25)
% 
xlabel('Time (fs)','Interpreter','tex');
ylabel('Kinetic Energy (Ha/atom)','Interpreter','tex');
box on; 
axis square;
set(gca,'LineWidth',lw);
set(gca,'FontSize',fsz);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');

timetot = cputime - tstart;
fprintf('*************************************************************** \n');
fprintf('\n');
fprintf(' Time to complete the FPMD Simulation = %20.3e s\n', timetot);
fprintf('\n');
fprintf('*************************************************************** \n');

