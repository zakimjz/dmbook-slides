
#*********************************************************************************************

#directory where the generated data will be stored
#data_directory = "D:/Gabi/Implementation/DataSets/KDD08/Uniform/Equal/Parallel/"
data_directory = "/tmp/"

#file with the generated data
file_output <- paste(data_directory, "Set1.txt", sep="")
file.create(file_output)

#file with the relevant attributes for each cluster
file_reldim <- paste(data_directory, "Set1_RelDim.txt", sep="")
file.create(file_reldim)

#file with the data in the format required by Biosphere
file_forBiosphere <- paste(data_directory, "Set1_forBiosphere.txt", sep="")
file.create(file_forBiosphere)

rfactor = 1
#statistical significance level
alpha_zero = 1.0E-5

#range of each attribute for the generated data
range_min = 0
range_max = 1

#number clusters
k = 5
#number attributes
d = 50
#number objects
n = 300*rfactor

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

# relevant attributes for each cluster
R <- matrix(rep(0,k*d), nrow = k, ncol = d)

#for each cluster, for each relevant attribute, we keep the defining range
Range_Min <- matrix(rep(0,k*d), nrow = k, ncol = d)
#for each cluster, for each relevant attribute, we keep the defining range
Range_Max <- matrix(rep(0,k*d), nrow = k, ncol = d)

###################################################################

# generate cluster 1
IDCluster = 1
Volume = 1
#cluster 1 on attribute 1 spans the interval [0.2,0.5]
Range_Min[IDCluster,1] = 0.2 
Range_Max[IDCluster,1] = 0.5
Volume = Volume * (Range_Max[IDCluster,1] - Range_Min[IDCluster,1])
#cluster 1 on attribute 2 spans the interval [0.1,0.4]
Range_Min[IDCluster,2] = 0.1
Range_Max[IDCluster,2] = 0.4
Volume = Volume * (Range_Max[IDCluster,2] - Range_Min[IDCluster,2])

# determine the critical value
ExpectedSupport = n * Volume
alpha_adjusted = alpha_zero / choose(d, length(rel_dim1))
theta_alpha = 0
bSatisfied = 0
sum = 0
while(bSatisfied == 0){

	#sum = sum + dpois(theta_alpha, ExpectedSupport)
	sum = sum + dbinom(theta_alpha, n, Volume)
	if (sum >= 1 - alpha_adjusted){
		bSatisfied = 1	
	}
	else{
		theta_alpha = theta_alpha + 1
	}	
}
print(theta_alpha)

# number points in cluster
#n1 = theta_alpha + 1
n1 = 60*rfactor
C1 <- matrix(nrow = n1, ncol = d)
for(j in 1:d){
	if (length(rel_dim1[rel_dim1 == j]) == 0){
		#irrelevant attribute
		C1[,j] <- runif(n1,min = range_min, max = range_max)
	}
	else{
		C1[,j] <- runif(n1, min = Range_Min[IDCluster,j], max = Range_Max[IDCluster,j])
		R[IDCluster,j] <- 1
	}
}
true_labels <- rep(1,n1)

###################################################################

# generate cluster 2
IDCluster <- 2
Volume = 1
Range_Min[IDCluster,3] = 0.2 
Range_Max[IDCluster,3] = 0.4
Volume = Volume * (Range_Max[IDCluster,3] - Range_Min[IDCluster,3])
Range_Min[IDCluster,4] = 0.1
Range_Max[IDCluster,4] = 0.4
Volume = Volume * (Range_Max[IDCluster,4] - Range_Min[IDCluster,4])

# determine the critical value
ExpectedSupport = n * Volume
alpha_adjusted = alpha_zero / choose(d, length(rel_dim2))
theta_alpha = 0
bSatisfied = 0
sum = 0
while(bSatisfied == 0){

	#sum = sum + dpois(theta_alpha, ExpectedSupport)
	sum = sum + dbinom(theta_alpha, n, Volume)
	if (sum >= 1 - alpha_adjusted){
		bSatisfied = 1	
	}
	else{
		theta_alpha = theta_alpha + 1
	}	
}
print(theta_alpha)

# number points in cluster
#n2 = theta_alpha + 1
n2 = 50*rfactor
C2 <- matrix(nrow = n2, ncol = d)
for(j in 1:d){
	if (length(rel_dim2[rel_dim2 == j]) == 0){
		#irrelevant attribute
		C2[,j] <- runif(n2,min = range_min, max = range_max)
	}
	else{
		C2[,j] <- runif(n2, min = Range_Min[IDCluster,j], max = Range_Max[IDCluster,j])
		R[IDCluster,j] <- 1
	}
}
true_labels <- c(true_labels, rep(2,n2))

###################################################################

# generate cluster 3
IDCluster <- 3
Volume = 1
Range_Min[IDCluster,5] = 0.2 
Range_Max[IDCluster,5] = 0.4
Volume = Volume * (Range_Max[IDCluster,5] - Range_Min[IDCluster,5])
Range_Min[IDCluster,6] = 0.1
Range_Max[IDCluster,6] = 0.3
Volume = Volume * (Range_Max[IDCluster,6] - Range_Min[IDCluster,6])

# determine the critical value
ExpectedSupport = n * Volume
alpha_adjusted = alpha_zero / choose(d, length(rel_dim3))
theta_alpha = 0
bSatisfied = 0
sum = 0
while(bSatisfied == 0){

	#sum = sum + dpois(theta_alpha, ExpectedSupport)
	sum = sum + dbinom(theta_alpha, n, Volume)
	if (sum >= 1 - alpha_adjusted){
		bSatisfied = 1	
	}
	else{
		theta_alpha = theta_alpha + 1
	}	
}
print(theta_alpha)

# number points in cluster
#n3 = theta_alpha + 1
n3 = 40*rfactor
C3 <- matrix(nrow = n3, ncol = d)
for(j in 1:d){
	if (length(rel_dim3[rel_dim3 == j]) == 0){
		#irrelevant attribute
		C3[,j] <- runif(n3,min = range_min, max = range_max)
	}
	else{
		C3[,j] <- runif(n3, min = Range_Min[IDCluster,j], max = Range_Max[IDCluster,j])
		R[IDCluster,j] <- 1
	}
}
true_labels <- c(true_labels, rep(3,n3))

###################################################################

# generate cluster 4
IDCluster <- 4
Volume = 1
Range_Min[IDCluster,7] = 0.2 
Range_Max[IDCluster,7] = 0.3
Volume = Volume * (Range_Max[IDCluster,7] - Range_Min[IDCluster,7])
Range_Min[IDCluster,8] = 0.1
Range_Max[IDCluster,8] = 0.4
Volume = Volume * (Range_Max[IDCluster,8] - Range_Min[IDCluster,8])

# determine the critical value
ExpectedSupport = n * Volume
alpha_adjusted = alpha_zero / choose(d, length(rel_dim4))
theta_alpha = 0
bSatisfied = 0
sum = 0
while(bSatisfied == 0){

	#sum = sum + dpois(theta_alpha, ExpectedSupport)
	sum = sum + dbinom(theta_alpha, n, Volume)
	if (sum >= 1 - alpha_adjusted){
		bSatisfied = 1	
	}
	else{
		theta_alpha = theta_alpha + 1
	}	
}
print(theta_alpha)

# number points in cluster
#n4 = theta_alpha + 1
n4 = 40*rfactor
C4 <- matrix(nrow = n4, ncol = d)
for(j in 1:d){
	if (length(rel_dim4[rel_dim4 == j]) == 0){
		#irrelevant attribute
		C4[,j] <- runif(n4,min = range_min, max = range_max)
	}
	else{
		C4[,j] <- runif(n4, min = Range_Min[IDCluster,j], max = Range_Max[IDCluster,j])
		R[IDCluster,j] <- 1
	}
}
true_labels <- c(true_labels, rep(4,n4))

###################################################################

# generate cluster 5
IDCluster <- 5
Volume = 1
Range_Min[IDCluster,9] = 0.2 
Range_Max[IDCluster,9] = 0.4
Volume = Volume * (Range_Max[IDCluster,9] - Range_Min[IDCluster,9])
Range_Min[IDCluster,10] = 0.3
Range_Max[IDCluster,10] = 0.6
Volume = Volume * (Range_Max[IDCluster,10] - Range_Min[IDCluster,10])

# determine the critical value
ExpectedSupport = n * Volume
alpha_adjusted = alpha_zero / choose(d, length(rel_dim5))
theta_alpha = 0
bSatisfied = 0
sum = 0
while(bSatisfied == 0){

	#sum = sum + dpois(theta_alpha, ExpectedSupport)
	sum = sum + dbinom(theta_alpha, n, Volume)
	if (sum >= 1 - alpha_adjusted){
		bSatisfied = 1	
	}
	else{
		theta_alpha = theta_alpha + 1
	}	
}
print(theta_alpha)

# number points in cluster
#n5 = theta_alpha + 1
n5 = 50*rfactor
C5 <- matrix(nrow = n5, ncol = d)
for(j in 1:d){
	if (length(rel_dim5[rel_dim5 == j]) == 0){
		#irrelevant attribute
		C5[,j] <- runif(n5,min = range_min, max = range_max)
	}
	else{
		C5[,j] <- runif(n5, min = Range_Min[IDCluster,j], max = Range_Max[IDCluster,j])
		R[IDCluster,j] <- 1
	}
}
true_labels <- c(true_labels, rep(5,n5))

###################################################################

#the noise points
#no_noise_points = n - (n1 + n2 + n3 + n4 + n5)
no_noise_points = 60*rfactor
no_noise_points = 0*rfactor
true_labels <- c(true_labels,rep(k+1,no_noise_points))

#generate the noise
N <- matrix(nrow = no_noise_points, ncol= d)
if (no_noise_points > 0){
	for(i in 1:no_noise_points){
		N[i,] <- runif(d, min = range_min, max = range_max)
	}
}

#generate the whole data set
X <- rbind(C1,C2,C3,C4,C5,N)


#******************************************************************************************************

# noise objects must be relabeled if they fall within the hyper-rectangle of a cluster

for(l in 1:k){ #for each cluster
	x <- 1 # keep in x the objects that fall within the range for each relevant attribute
	for(j in 1:d){ 	
		if (R[l,j] == 1){ #for each relevant attribute j
			x <- x & (X[,j] >= Range_Min[l,j]) & (X[,j] <= Range_Max[l,j]) 			
		}
	}
	s <- sort(x, index.return = TRUE)
	new_cluster_members <- s$ix[s$x == 1]
	for(j in 1:length(new_cluster_members)){
		if (true_labels[new_cluster_members[j]] == (k+1)){ #noise point
			true_labels[new_cluster_members[j]] = l
		}
	}
}


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


