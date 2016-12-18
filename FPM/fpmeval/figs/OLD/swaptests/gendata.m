
function D = gendata(n,d,c,e); 

%% generate 0/1 clustered data
%% n points
%% d dimensions
%% c clusters
%% e noise around each cluster

ppc = round(n/c); 
D = []; 
for i=1:c
  ci = round(rand(1,d)); 
  D = [ D ; ci(ones(ppc,1),:) ]; 
end
D = mod(randgd(size(D,1),size(D,2),e) + D, 2); 

