function Result = count3d(X,v1l,v1u,v2l,v2u,v3l,v3u)

% X is the 3d matrix,
% v1 and v2 is the range

c=0;
for i=1:60
  if X(i,1) < v1u && X(i,1) >= v1l && X(i,2) < v2u && X(i,2) >= v2l ...
        && X(i,3) < v3u && X(i,3) >= v3l
   c= c+1;
  end
end
Result=c;
