
function A = randgd(m,n,d); 

%% random m x n binary matrix with given density d

A = rand(m,n); 
[i,j] = find(A<=d); 
A = full(spconvert([ i j ones(size(i)) ; m n 0 ])); 
