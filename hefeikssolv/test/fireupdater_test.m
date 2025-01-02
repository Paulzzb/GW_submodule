function [v, dt, alpha, cut] = fireupdater_test(it, dt, F, v, cut, alpha, pars)
% 
% Subhajit Banerjee
% June 2017
% SSG@CRD, LBNL
% 
if (nargin == 6)
    pars.fDec = 0.5;
    pars.fInc = 1.1;
    pars.nMin = 5;
    pars.alphaStart = 0.1;
    pars.fAlpha = 0.99;
end

% 
fDec            = pars.fDec;
fInc            = pars.fInc;
nMin            = pars.nMin;
alphaStart      = pars.alphaStart;
fAlpha          = pars.fAlpha;
% 
dtMax = 10.0*dt;
% 
P = dot(F, v);      % Power
hatF = F/norm(F);   % Fhat (the unit vector)
%
% FIRE velocity update formula
%
v = (1.0 - alpha)*v + alpha*hatF*norm(v);
disp('alpha input = ')
disp(alpha);

disp('dt input = ')
disp(dt);

if P < 0.0
    v = 0.0;                    % Reset velosity to 0
    cut = it;                   % cut <-- iter # (gets updated everytime P <= 0)
    dt = dt*fDec;               % decrease dt
    alpha = alphaStart;         % reset alpha to alpha_start
elseif (P >= 0.0) && ((it - cut) > nMin)
    dt = min(dt*fInc, dtMax);   % Update/INCREASE dt
    alpha = fAlpha*alpha;       % update alpha
end
disp('it = ')
disp(it)
disp('cut = ')
disp(cut)
disp('alpha output = ')
disp(alpha);
disp('dt output = ')
disp(dt);

return
end