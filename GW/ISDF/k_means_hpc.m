function ind_mu = k_means(rk, weight, options)
% Inputs:
%	  rk: number of categories in K-means
%	  weight: weight function
%	  points(in varargin): a set of points {ri} i=1,2...Nr
%
% Outputs:
%	  ind_mu: store indice of centroids data

option=options.isdfoptions;
n1=option.sys.n1;n2=option.sys.n2;n3=option.sys.n3;
[I,J,K] = ndgrid((0:n1-1)/n1-((0:n1-1) >= n1/2), ...
	(0:n2-1)/n2-((0:n2-1) >= n2/2), ...
	(0:n3-1)/n3-((0:n3-1) >= n3/2));
points=reshape(cat(4,I,J,K),[],3)*option.sys.supercell;
atompos=option.sys.xyzlist;

% Nr: number of points; nd:dimension of points
[Nr, nd] = size(points);
assert(rk<=Nr,'The number of points is less than rk.')

rng(option.seed,'twister');
cluster = zeros(Nr, 1);
ind_mu = zeros(rk, 1);
dist = zeros(Nr, rk);

% initial centroids index
tic
switch option.init
	case 'sequence'
		lastCentroids_ind_mu=round(1:floor(Nr/rk):Nr);
		lastCentroids_ind_mu=lastCentroids_ind_mu(1:rk)';
	case 'random'
		lastCentroids_ind_mu = randperm(Nr, rk)';
	case 'wrs'
		% Weighted random sampling
		weight = abs(weight);
		lastCentroids_ind_mu = datasample(1:Nr,rk,'Replace',false,'Weights',weight);
	case 'kmeans++_false'
		lastCentroids_ind_mu(1)=randi([1,Nr],1);
		for i=2:rk
			dist(:,i-1)=sum((points-points(lastCentroids_ind_mu(i-1),:)).^2,2);
			[min_dist, ~] = min(dist(:,1:i-1), [], 2);
			weighted_dist=min_dist.*weight;
			lastCentroids_ind_mu(i)=randsample(1:Nr,1,true,weighted_dist);
		end
	case 'kmeans++'
		lastCentroids_ind_mu(1)=randi([1,Nr],1);
		for i=2:rk
			dist(:,i-1)=sum((points-points(lastCentroids_ind_mu(i-1),:)).^2,2);
			[min_dist, ~] = min(dist(:,1:i-1), [], 2);
			weighted_dist=min_dist.*weight;
			[sorted_dist,I] = sort(weighted_dist,'descend');
			rand_dist=0.9*rand*sum(sorted_dist);
			sum_dist=0;
			for j=1:Nr-1
				sum_dist=sum_dist+sorted_dist(j);
				if sum_dist>=rand_dist
					lastCentroids_ind_mu(i)=I(j);
					break
				end
			end
		end
	case 'atom'
		for i=1:size(atompos,1)
			dist(:,i)=sum((points-atompos(i,:)).^2,2);
			[~,I]=min(dist(:,i));
			lastCentroids_ind_mu(i)=I;
		end
		for i=size(atompos,1)+1:rk
			dist(:,i-1)=sum((points-points(lastCentroids_ind_mu(i-1),:)).^2,2);
			[min_dist, ~] = min(dist(:,1:i-1), [], 2);
			weighted_dist=min_dist.*weight;
			[sorted_dist,I] = sort(weighted_dist,'descend');
			rand_dist=0.9*rand*sum(sorted_dist);
			sum_dist=0;
			for j=1:Nr-1
				sum_dist=sum_dist+sorted_dist(j);
				if sum_dist>=rand_dist
					lastCentroids_ind_mu(i)=I(j);
					break
				end
			end
		end
	otherwise
		error('Unknown initialization method: %s', option.init);
end
%fprintf('kmeans init time: %f\n',toc);
% the new centroids index
newCentroids_ind_mu = zeros(rk, 1);

max_iteration=100;
iteration = 1;
%distance between grids and centroids
tic
while true
	if iteration == 1
		timefirst = tic;
	end
	for nk=1:rk
		dist(:,nk) = sum((points-points(lastCentroids_ind_mu(nk),:)).^2,2);
	end
	[~, index_min] = min(dist, [], 2);
	shifted_points=points;
	 
	% calulate centroids
	parfor mk=1:rk
		%grids belong to the same cluster
		cluster=find(index_min==mk);
		if isempty(cluster)
			if newCentroids_ind_mu(mk) == 0
				newCentroids_ind_mu(mk) = 1;
			end
		  continue;
		end
		total_pointsMutilweight = sum(shifted_points(cluster, :).*weight(cluster),1);
		total_weight = sum(weight(cluster),1);
		Centroids = total_pointsMutilweight/ total_weight;
		
		% find the nearest point as estimated centroids
		distPoint2Centroid = sum((Centroids-shifted_points(cluster,:)).^2, 2);
		[~, I] = min(distPoint2Centroid);
		newCentroids_ind_mu(mk,1) = cluster(I);
	end
	
	%The calculation will be stopped if one of	two conditions is satisfied
	%reach the number of scheduled max_iteration or  the centroids change within convergence criteria
	if max_iteration ~= 0 && iteration == max_iteration
		warning('K-Means reach max iterations!');
		ind_mu = newCentroids_ind_mu;
		break;
	elseif	norm(newCentroids_ind_mu-lastCentroids_ind_mu) < 0.1 %Convergence criteria
		ind_mu = newCentroids_ind_mu;
		%fprintf('kmeans iteration time: %f\n',toc);
		break;
	end
	if iteration == 1
		timefirst = toc(timefirst);
	  fprintf("Time for first iteration = %6.2f.\n", timefirst);
		fprintf("Max Time = %6.2f.\n", timefirst * max_iteration);
	end
	iteration = iteration + 1;
	lastCentroids_ind_mu = newCentroids_ind_mu;
end
end
