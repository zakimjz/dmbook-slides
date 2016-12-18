function Result = count1d(X,a,v1l,v1u)

% X is the 3d matrix, a and b are the cols for 2d proj
% v1 and v2 is the range

c=0;
for i=1:60
  if X(i,a) < v1u && X(i,a) >= v1l
   c= c+1;
  end
end
Result=c;
