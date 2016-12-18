X <- read.table('2d-3c-no7.dat')
Y <- read.table('2d-3c-no7-partition.out')
X2 <- read.table('2d-3c-no19.dat')
Y2 <- read.table('2d-3c-no19-partition.out')


D <- cbind(X,Y)
DD <- cbind(X2,Y2)
D1 <- subset(D, D[,4]==1)
D2 <- subset(D, D[,5]==1)
D3 <- subset(D, D[,6]==1)
DD1 <- subset(DD, DD[,3]==1)
DD2 <- subset(DD, DD[,4]==1)
DD3 <- subset(DD, DD[,5]==1)


print(c(nrow(D), nrow(D1), nrow(D2), nrow(D3)))
print(c(nrow(DD), nrow(DD1), nrow(DD2), nrow(DD3)))
write.table(D1[,1:3], 'n7-2d-d1.dat', row.names=F, col.names=F)
write.table(D2[,1:3], 'n7-2d-d2.dat', row.names=F, col.names=F)
write.table(D3[,1:3], 'n7-2d-d3.dat', row.names=F, col.names=F)
write.table(DD1[,1:2], 'n19-2d-d1.dat', row.names=F, col.names=F)
write.table(DD2[,1:2], 'n19-2d-d2.dat', row.names=F, col.names=F)
write.table(DD3[,1:2], 'n19-2d-d3.dat', row.names=F, col.names=F)
