
function [orig,permuted,sw] = clusteringtest(Dinp,k,nosamples,L);

%% Clustering test 

%% Input 
%%   Dinp: input 0-1 matrix
%%   k: number of clusters to be tested
%%   nosamples: number of samples (swapped datasets) to be drawn
%%   L: number of swaps to be attempted

%% Output
%%   orig: error on the original data
%%   permuted: vector of errors on permuted data
%%   sw: vector of actual swaps took place

D=Dinp;

if (nargin<=2)
  nosamples = 100; 
end
  
sw = zeros(1,nosamples);

if (size(D,1)<=k)
  err = 0; 
else
  [centroids,clusters,err] = kmeans(D,k,'EmptyAction','drop');  
end
orig = sum(err);

for i=1:nosamples
  if (size(D,1)<=k)
    err = 0; 
    swaps = 0; 
  else
    if (nargin>=4)
      [D,swaps] = swap(D,L);
    else 
      [D,swaps] = swap(D);
    end
    [centroids,clusters,err] = kmeans(D,k,'EmptyAction','drop');
  end
  permuted(i) = sum(err);
  sw(i) = swaps;
end




