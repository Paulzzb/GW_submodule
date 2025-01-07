% This code is used to test parfor in kmeans main-loop.
% The main loop contains
n1 = 20; n2 = 20; n3 = 20;
nx = n1*n2*n3; nc = n1*n2;

[I,J,K] = ndgrid((0:n1-1)/n1-((0:n1-1) >= n1/2), ...
  (0:n2-1)/n2-((0:n2-1) >= n2/2), ...
  (0:n3-1)/n3-((0:n3-1) >= n3/2));
points=reshape(cat(4,I,J,K),[],3);

xlist = points; clist = zeros(nc, 3);
cinds = zeros(nc, 1); xclass = zeros(nx, 1);

cinds = randi(nx, nc, 1);
clist = xlist(cinds,:);

flagconv = 1;
while flagconv
  % 1. Classify :
  for i = 1:nx
    x2c = sum((xlist(i, :) - clist(:, :)).^2, 2);
    [~, ind] = min(x2c);
    xclass(i) = ind;
  end
  % 2. Find center: 
  for i = 1:nc
    xindinc = find(xclass == c);
    nxinc = length(xindinc);
    xinc = xlist(xindinc, :);
    weighx = weight(xindinc);
    c = sum(xinc .* weighx, 2) / weightx;
    
    % Calculate ||x_k - c_k||, find the smallest one.
    distxc = sum((xlist(xindinc(j)) - c).^2, 2);
    [~, ind] = min(distxc);
    cinds(i) = ind
  end
  
  clist = xlist(cinds, :)
end % while
     
