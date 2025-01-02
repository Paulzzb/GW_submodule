%
% construct Si8 
%
kssolvpptype('pz-hgh', 'UPF');
%
% 1. construct atoms
%
a1 = Atom('Si');
atomlist = [a1, a1, a1, a1, a1, a1, a1, a1];
%
% 2. set up supercell
%
C = 10.216*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
    0.000000000         0.000000000         0.000000000
    0.000000000         0.500000000         0.500000000
    0.500000000         0.000000000         0.500000000
    0.500000000         0.500000000         0.000000000
    0.750000000         0.250000000         0.750000000
    0.250000000         0.250000000         0.250000000
    0.250000000         0.750000000         0.750000000
    0.750000000         0.750000000         0.250000000
];
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule('supercell',C, 'atomlist',atomlist, 'xyzlist' ,xyzlist, ...
    'ecut',12.5, 'name','Si8' );

% 
% Set up MD run
% 
nAtoms = sum(mol.natoms);
ksopts = setksopt();
ksopts.maxscfiter = 200;
ksopts.betamix = 0.08;
ksopts.mixdim = 6;
ksopts.scftol = 1e-6;

ksopts = setksopt(ksopts,'mixtype','broyden','eigmethod','eigs');
[mol,H0,X0,info] = scf(mol, ksopts);

nMdSteps = 1000;

[X, V, F, eTot, eKin, eFree, entropy, INFO] = md(mol, H0, X0, nMdSteps, ksopts);

dlmcell('Si8_MD_runSummary.txt', INFO);

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
plot(1:nMdSteps, eKin/nAtoms, 'b', 'LineWidth', 1.25)
% 
xlabel('Time (fs)','Interpreter','tex');
ylabel('Kinetic Energy (Ha/atom)','Interpreter','tex');
box on; 
axis square;
set(gca,'LineWidth',lw);
set(gca,'FontSize',fsz);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');
% 
% Plot Free energy:
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
plot(1:nMdSteps, eFree/nAtoms, 'k', 'LineWidth', 1.25);
% 
xlabel('Time (fs)','Interpreter','tex');
ylabel('Free Energy (Ha/atom)','Interpreter','tex');
box on; 
axis square;
set(gca,'LineWidth',lw);
set(gca,'FontSize',fsz);
set(gca,'FontWeight','Bold');
set(gcf,'color','w');
% 
% Plot Entropy:
% 
% figure
% set(gcf,'InvertHardcopy','on');
% set(gcf,'PaperUnits', 'inches');
% pos = get(gcf, 'Position');
% set(gcf, 'Position', [pos(1) pos(2) width*1000, height*1000]);
% set(gca, 'FontSize', fsz, 'LineWidth', alw);
% set(get(gca,'xlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% set(get(gca,'ylabel'),'FontSize', fsz, 'FontWeight', 'Bold');
% set(get(gca,'title'),'FontSize', fsz, 'FontWeight', 'Bold');
% hold on;
% % 
% plot(1:nMdSteps, entropy/nAtoms, 'm', 'LineWidth', 1.25);
% % 
% xlabel('Time (fs)','Interpreter','tex');
% ylabel('Entropy (Ha/K/atom)','Interpreter','latex');
% box on; 
% axis square;
% set(gca,'LineWidth',lw);
% set(gca,'FontSize',fsz);
% set(gca,'FontWeight','Bold');
% set(gcf,'color','w');