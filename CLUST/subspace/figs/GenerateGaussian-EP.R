
#*********************************************************************************************

#directory where the generated data will be stored
#data_directory = "D:/Gabi/Implementation/DataSets/KDD08/Gaussian/Equal/Parallel/"
data_directory = "~/Implementation/DataSets/KDD08/Gaussian/Equal/Parallel/"

#file with the generated data
file_output <- paste(data_directory, "Set1.txt", sep="")
file.create(file_output)

#file with the relevant attributes for each cluster
file_reldim <- paste(data_directory, "Set1_RelDim.txt", sep="")
file.create(file_reldim)

#file with the data in the format required by Biosphere
file_forBiosphere <- paste(data_directory, "Set1_forBiosphere.txt", sep="")
file.create(file_forBiosphere)

#statistical significance level
#alpha_zero = 1.0E-5

#interval standard deviation
var_min = 0.004
var_max = 0.008

#range of each attribute for the generated data
range_min = 0
range_max = 1

#number clusters
k = 5
#number attributes
d = 50
#number objects
n = 300

#cluster 1
#relevant attributes
rel_dim1 <- c(1,2)

#cluster 2
#relevant attributes
rel_dim2<- c(3,4)

#cluster 3
#relevant attributes
rel_dim3<- c(5,6)

#cluster 4
#relevant attributes
rel_dim4<- c(7,8)

#cluster 5
#relevant attributes
rel_dim5<- c(9,10)

#******************************************************************************************************

#write a specific header in the file for Biosphere
cat("HEADER",file = file_forBiosphere, sep = "", append = TRUE)
cat("\t",file = file_forBiosphere, sep = "", append = TRUE)
for(j in 1:d){
	cat("N",file = file_forBiosphere, sep = "", append = TRUE)
	cat("\t",file = file_forBiosphere, sep = "", append = TRUE)
}
cat("CLASS",file = file_forBiosphere, sep = "", append = TRUE)
cat("\n",file = file_forBiosphere, sep = "", append = TRUE)

#write relevant attributes 
write_to_file <- function(id_cluster, rel_dim , file_reldim, file_forBiosphere ){

	cat("Cluster",file = file_reldim, sep = " ", append = TRUE)
	cat(id_cluster,file = file_reldim, sep = " ", append = TRUE)
	#cat(" ",file = file_output, sep = "", append = TRUE)
	
	cat("c",file = file_forBiosphere, sep = "", append = TRUE)
	cat(id_cluster-1,file = file_forBiosphere, sep = "", append = TRUE)
	cat("\t",file = file_forBiosphere, sep = "", append = TRUE)

	for(j in 1:d){
		if (length(rel_dim[rel_dim==j])==0){
			#j is not a relevant dimension
			cat(" ",file = file_reldim, sep = "", append = TRUE)
			cat(0,file = file_reldim, sep = " ", append = TRUE)
			cat(0,file = file_forBiosphere, sep = "", append = TRUE)
			cat("\t",file = file_forBiosphere, sep = "", append = TRUE)
		}
		else{
			cat(" ",file = file_reldim, sep = "", append = TRUE)
			cat(1,file = file_reldim, sep = " ", append = TRUE)
			cat(1,file = file_forBiosphere, sep = "", append = TRUE)
			cat("\t",file = file_forBiosphere, sep = "", append = TRUE)
		}
	}
	cat("\n",file = file_reldim, sep = "", append = TRUE)
	cat("\n",file = file_forBiosphere, sep = "", append = TRUE)		
}

write_to_file(1,rel_dim1,file_reldim,file_forBiosphere)
write_to_file(2,rel_dim2,file_reldim,file_forBiosphere)
write_to_file(3,rel_dim3,file_reldim,file_forBiosphere)
write_to_file(4,rel_dim4,file_reldim,file_forBiosphere)
write_to_file(5,rel_dim5,file_reldim,file_forBiosphere)


#******************************************************************************************************

library(MASS)
generate_cluster <- function(n,rel_dim){

	#fill the relevant dimensions
	no_rel_dim <- length(rel_dim)
	#anchor point
	mu <- runif(no_rel_dim, min = range_min, max = range_max)
	#covariance matrix
	sigma <- diag(,nrow = no_rel_dim, ncol = no_rel_dim)
	for(i in 1:no_rel_dim){
		var_i <- runif(1, var_min, var_max)
		sigma[i,i] <- var_i
	}

	R <- mvrnorm(n,mu,sigma,empirical=TRUE)
	
	#translate to the given ranges
	for(i in 1:no_rel_dim){
		if (min(R[,i])<range_min){
			R[,i] <- R[,i] - min(R[,i])
		}
		else{
			if (max(R[,i])>range_max){
				R[,i] <- R[,i] - (max(R[,i])-range_max)
			}
			
		}
	}
	

	#the points
	P <- matrix(nrow = n, ncol = d)
	contor <- 0
	for(i in 1:d){
		if (length(rel_dim[rel_dim==i])==0){
			#irrelevant dimension
			u <- runif(n,min = range_min, max = range_max)
			P[,i] <- u
		}
		else{
			contor = contor + 1
			P[,i] <- R[,contor]		
		}
	}
	P
}

###################################################################

# generate cluster 1
n1 = 60
C1 <- generate_cluster(n1,rel_dim1)
true_labels <- rep(1,n1)

###################################################################

# generate cluster 2
n2 = 50
C2 <- generate_cluster(n2,rel_dim2)
true_labels <- c(true_labels, rep(2,n2))

###################################################################

# generate cluster 3
n3 = 40
C3 <- generate_cluster(n3,rel_dim3)
true_labels <- c(true_labels, rep(3,n3))

###################################################################

# generate cluster 4
n4 = 40
C4 <- generate_cluster(n4,rel_dim4)
true_labels <- c(true_labels, rep(4,n4))

###################################################################

# generate cluster 5
n5 = 50
C5 <- generate_cluster(n5,rel_dim5)
true_labels <- c(true_labels, rep(5,n5))

###################################################################

#the noise points
no_noise_points = n - (n1 + n2 + n3 + n4 + n5)
true_labels <- c(true_labels,rep(k+1,no_noise_points))

#generate the noise
N <- matrix(nrow = no_noise_points, ncol= d)
for(i in 1:no_noise_points){
	N[i,] <- runif(d, min = range_min, max = range_max)
}

#generate the whole data set
X <- rbind(C1,C2,C3,C4,C5,N)


#******************************************************************************************************


#*****************************************************************************************************

## min-max normalization on each attribute
n <- length(X[,1])
for(j in 1:d){
	MaxVal <- max(X[,j])
	MinVal <- min(X[,j])
	for(i in 1:n){
		if (MaxVal != MinVal){
			X[i,j] <- (X[i,j] - MinVal)/(MaxVal - MinVal)
		}
		else{
			X[i,j] <- 0
		}
	}
	
}

#******************************************************************************************************

# write the data matrix to the files
#generate random ordering 
n <- length(X[,1])
order <- sample(1:n)
#order <- c(1:n)
for (i in 1:n){
	#file output
	cat(X[order[i],],file = file_output, sep = " ", append = TRUE)
	cat(" ",file = file_output, sep = "", append = TRUE)
	cat(true_labels[order[i]],file = file_output, sep = "", append = TRUE)
	cat("\n",file = file_output, sep = "", append = TRUE)

	#file for Biosphere
	cat("\t",file = file_forBiosphere, sep = "", append = TRUE)
	cat(X[order[i],],file = file_forBiosphere, sep = "\t", append = TRUE)
	cat("\t",file = file_forBiosphere, sep = "", append = TRUE)
	cat("c",file = file_forBiosphere, sep = "", append = TRUE)
	cat(true_labels[order[i]]-1,file = file_forBiosphere, sep = "", append = TRUE)
	cat("\n",file = file_forBiosphere, sep = "", append = TRUE)
}


