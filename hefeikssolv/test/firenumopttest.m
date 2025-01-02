% function [xValRec, fValRec, forceValRec] = firenumopttest(plotFlag, objFuncHn)
function firenumopttest(plotFlag, objFuncHn)
clc
% 
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL
% 
if nargin == 0
    plotFlag = 'yes';
    objFuncHn = @pes;
    % objFuncHn = @quadfunc;
%     objFuncHn = @Rosenbrockfunc;
end
%
% Source:
% https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.97.170201
%
% Qualitatively they should satisfy: nMin larger than 1 (at
% least a few smooth steps after freezing); fInc larger than but
% near to one (avoid too fast acceleration); fDec smaller than
% 1 but much larger than zero (avoid too heavy slowing down),
% alphaStart larger than, but near to zero (avoid too heavy damping);
% fAlpha smaller than, but near to one (mixing is efficient some time
% after restart).
%
pars.fDec = 0.5;
pars.fInc = 1.1;
pars.nMin = 3;
pars.alphaStart = 0.1;
pars.fAlpha = 0.99;

dt = 1/100;
maxit = 1000;
TOL = 1.0E-6;
mass = 1.0;
alpha = pars.alphaStart;
% 
% Initialization
% 
% x0 = 0.01*ones(2,1);
% x0 = zeros(2,1) + 1.0E-06;
% x0 = [0.0000001078; 1.2283902077] + 0.001;    % matlab minimizer
x0 = rand(2,1);
v0 = zeros(size(x0));                  % velocity initialized to 0 for FIRE
xIni = x0;

fRec = zeros(maxit,1);
xRec = cell(maxit,1);
forceValRec = cell(maxit,1);
%  
[E0, fValPrime] = objFuncHn(x0);
F0 = -fValPrime;
% 
% First velocity update with initial velocity = 0
%
[x, v, E, F] = velverletgen(x0, v0, F0, dt, mass, objFuncHn);

it = 1;
fRec(it) = E0;
xRec{it} = x0;
forceValRec{it} = F0;

cut = 0;

while (it <= maxit)
    
    if it == maxit
        fprintf(' FIRE did not converge! \n');
        fprintf(' Reached maximum number ( %-5d) of iterations! \n', it);
        break
    end
    
    % Check for convergence
    checkConv = norm(F);
%     checkConv = max(abs(F));

    if checkConv < TOL
        break
    else
        % Update v, dt, alpha using FIRE
        [v, dt, alpha, cut] = fireupdater_test(it, dt, F, v, cut, alpha, pars);
        
        %
        % Rerun MD
        %
        x0 = x;
        v0 = v;
        F0 = F;
        it = it + 1;
                
        xRec{it} = x;
        fRec(it) = E;
        forceValRec{it} = F;
        
        [x, v, E, F] = velverletgen(x0, v0, F0, dt, mass, objFuncHn);

    end
end

maxNnzEl = max(find(fRec~=0));
fValRec =  fRec(1:maxNnzEl);

xRecMat = [xRec{:}];
xValRec = xRecMat(:,1:maxNnzEl);

forceValRecMat = [forceValRec{:}];
forceValRec = forceValRecMat(:,1:maxNnzEl);

if (strcmpi(plotFlag, 'yes') == 1),
    
    xLow = min(xValRec(1,:));
    xUp = max(xValRec(1,:));
    
    yLow = min(xValRec(2,:));
    yUp = max(xValRec(2,:));
    
    if abs(xLow) < 1.0 && abs(xUp) < 1.0
        xLow = xLow*3.0;
        xUp = xUp*3.0;
    end
    if abs(yLow) < 1.0 && abs(yUp) < 1.0
        yLow = yLow*3;
        yUp = yUp*3;
    end
    
    % Resolution of the surface plot
    nPts = 100;
    
    xSampl = linspace(xLow, xUp, nPts);
    ySampl = linspace(yLow, yUp, nPts);
    
    [X, Y] = meshgrid(xSampl, ySampl);
    
    EVAL = zeros(size(X));
    FVAL = cell(size(X));

    for i = 1:size(X,1)
        for j = 1:size(Y,1),
            [tmp1, tmp2] = objFuncHn([X(i,j); Y(i,j)]);
            EVAL(i,j) = tmp1;
            FVAL{i,j} = -tmp2;
        end
    end
    
    close all
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
    set(get(gca,'zlabel'),'FontSize', fsz, 'FontWeight', 'Bold');
    set(get(gca,'title'),'FontSize', fsz, 'FontWeight', 'Bold');
    hold on;
    
    surf(X, Y, EVAL)
    
    axis square;
    colormap jet
    set(gca,'LineWidth',lw);
    set(gca,'FontSize',fsz);
    set(gca,'FontWeight','Bold');
    set(gcf,'color','w');
    set(gca,'CameraPosition',[10 -1 10]);
    view (-140,18);
    lighting phong
    shading interp    
    camlight left
    box on;
    hold on;
end

if (strcmpi(plotFlag, 'yes') == 1),
    plot3(xValRec(1,:), xValRec(2,:), fValRec, 'ro-');
end

% Options for BFGS
options = optimoptions(@fminunc,'Display','iter','Algorithm','quasi-newton');

% Options for trust-region
% optimoptions('fminunc','Display','iter','Algorithm','trust-region','SpecifyObjectiveGradient',true);  % Works for latest version of MATLAB
% options = optimoptions('fminunc','Display','iter','Algorithm','trust-region','GradObj','on');           % Works for MATLAB2013

if it == maxit
    [minF, ind] = min(fValRec);
    fprintf('\n');
    fprintf(' The infimum achieved is = %-16.10f \n', minF);
    fprintf(' The infimum is achieved at point (argmin) [%-16.10f, %-16.10f] \n',xValRec(1,ind), xValRec(2,ind));
    fprintf(' in %-5d iterations \n', ind); 
    fprintf('\n');
    [xValExact, fValExact, ~, info] = fminunc(objFuncHn, xIni, options);
    fprintf(' Optimization routine from MATLAB gives: %-16.10f \n', fValExact);
    fprintf(' This is achieved at point (argmin) [%-16.10f, %-16.10f] \n',xValExact(1), xValExact(2));
    fprintf(' Now Exiting \n')
    fprintf('\n');
    return
else
    fprintf('\n');
    fprintf(' FIRE converged in %-5d steps \n', it);
    fprintf(' The local minumum achieved is = %-16.10f \n', fValRec(end));
    fprintf(' This is achieved at point (argmin) [%-16.10f, %-16.10f] \n',xValRec(1,end), xValRec(2,end));
    fprintf('\n');
    
    fprintf('\n');
    [xValExact, fValExact, ~, info] = fminunc(objFuncHn, xIni, options);
    fprintf(' Optimization routine from MATLAB gives: %-16.10f \n', fValExact);
    fprintf(' This is achieved at point (argmin) [%-16.10f, %-16.10f] \n',xValExact(1), xValExact(2));
    fprintf(' Now Exiting \n')
    fprintf('\n');
    fprintf(' Now Exiting \n')
    fprintf('\n');
    return
end

end

function [fVal, fValPrime] = pes(x)
% 
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL
% 
argmt = pi*norm(x) + 0.5*atan(x(2)/x(1));
fVal = sin(argmt) + (norm(x))^2/10.0;
fValPrime = [cos(argmt)*(pi*x(1)/norm(x) - 0.5*x(2)/(norm(x))^2) + x(1)/5.0; 
             cos(argmt)*(pi*x(2)/norm(x) + 0.5*x(1)/(norm(x))^2) + x(2)/5.0];
end

function [f, fPrime] = Rosenbrockfunc(x, beta)
% 
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL
% 
if nargin == 1,
    beta = 10.0;
end

% Always input length(x) = Even number!
% The extended Rosenbrock function
n=length(x);
sum = 0;
gradf = zeros(n,1);
for i=1:n/2
    sum = sum + (beta*(x(2*i) - x(2*i-1).^2).^2 + (1 - x(2*i-1)).^2);
end
f = sum;
for j=2:2:n
    gradf(j) = gradf(j) + 2*beta*(x(j) - x(j-1).^2);
end
gradf(1) = gradf(1) - 4*beta*(x(2) - x(1).^2).*x(1) - 2*(1-x(1));
for k=3:2:(n-1)
    gradf(k) = gradf(k) - 4*beta*(x(k+1) - x(k).^2).*x(k) - 2*(1-x(k));
end
fPrime = gradf;
end

function [E, F] = quadfunc(x)
% 
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL
% 
E = dot(x, x);
F = 2.0*x;

end

function [x, v, fVal, F] = velverletgen(x0, v0, F0, dt, mass, objFuncHn)
% 
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL
% 
% INPUT:
% x0        = Initial Position vector
% F0        = Initial Force vector (negative gradient of the 
%                                   objective function)
% v0        = Initial velocity
% dt        = Time step
% objFuncHn = function handle pointing to function which evaluates the
%             objective function and its derivatives                                 

%
% Position update in Velocity Verlet
% 
x = x0 + v0*dt + 0.5*dt*dt*F0/mass;         % (F0/mass) is accln

% 
% Just one objective function evaluation:
% 
[fVal, fValPrime] = objFuncHn(x);      % [fVal, fValPrime] = objFuncHn(x);

F = -fValPrime;                        % force = -gradient

% Velocity update in Velocity Verlet
v = v0 + 0.5*(F0 + F)*dt/mass;

end




