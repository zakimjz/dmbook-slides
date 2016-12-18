X = load ('-ascii', 'subspace-db.txt');
C = [0.5, 1.5, 2.5, 3.5, 4.5]; %centers for binning

hist(X(:,1),C);
hist(X(:,2),C);
hist(X(:,3),C);

count1d(X,1,0,1)
count1d(X,1,1,2)
count1d(X,1,2,3)
count1d(X,1,3,4)
count1d(X,1,4,5)

count1d(X,2,0,1)
count1d(X,2,1,2)
count1d(X,2,2,3)
count1d(X,2,3,4)
count1d(X,2,4,5)

count1d(X,3,0,1)
count1d(X,3,1,2)
count1d(X,3,2,3)
count1d(X,3,3,4)
count1d(X,3,4,5)


count2d(X,1,2,0,1,0,1)
count2d(X,1,2,0,1,3,4)
count2d(X,1,2,0,1,4,5)
count2d(X,1,2,3,4,0,1)
count2d(X,1,2,3,4,3,4)
count2d(X,1,2,3,4,4,5)

count2d(X,1,3,0,1,0,1)
count2d(X,1,3,0,1,1,2)
count2d(X,1,3,0,1,3,4)
count2d(X,1,3,3,4,0,1)
count2d(X,1,3,3,4,1,2)
count2d(X,1,3,3,4,3,4)

count2d(X,2,3,0,1,0,1)
count2d(X,2,3,0,1,1,2)
count2d(X,2,3,0,1,3,4)
count2d(X,2,3,3,4,0,1)
count2d(X,2,3,3,4,1,2)
count2d(X,2,3,3,4,3,4)
count2d(X,2,3,4,5,0,1)
count2d(X,2,3,4,5,1,2)
count2d(X,2,3,4,5,3,4)

count3d(X,0,1,0,1,3,4)
count3d(X,0,1,4,5,1,2)
