## The original loop-based functions of lefse 0.5, kept verbatim as a
## reference implementation for the tests. Not exported.
lefse05 <- new.env()
local({
wt.sd <- function(x, wt) { s <- is.finite(x*wt); wt <- wt[s]; x <- x[s]; xbar <- sum(wt*x)/sum(wt); sqrt(sum(wt*(x-xbar)^2)*(sum(wt)/(sum(wt)^2-sum(wt^2)))) }
tips <- function(phy, node) geiger::tips(phy, node)
treedata <- function(...) geiger::treedata(...)
pd <- function(...) picante::pd(...)
Fmntd <-
function(dist.mat, my.sample){


Fmntd.sub = function(x){ 
		
## Get the names of the species present in a 
## community.
	com.names = names(x[x > 0])
		
## Make the community phylogenetic distance 
## matrix by extracting those rows and columns 
## that have species present in our community.
	my.com.dist = dist.mat[com.names, com.names]
		
## Set all diagonal values to NA so that the 
## zeros for conspecific comparisons do not 
## interfere with our calculation of nearest 
## neighbors.
	diag(my.com.dist) = NA
		
## Use apply() to calculate the minimum value in 
## each row of the community phylogenetic 
## distance matrix and take a mean of those 
## values.
	mean(apply(my.com.dist, MARGIN = 1, min, na.rm=T), na.rm=T)
	
}

apply(my.sample, MARGIN = 1, Fmntd.sub)

}

Fmntd.a <-
function(dist.mat, my.sample){

Fmntd.a.sub = function(x){ 
		
## Get the names of species present in the 
## community.
	com.names = names(x[x > 0])
		
## Make the community phylogenetic distance 
## matrix using the names of the species present ## in the community.
	my.com.dist = dist.mat[com.names, com.names]
		
## Place NA values in the diagonals
	diag(my.com.dist) = NA
		
## Calculate a mean of the minimum values in each 
## row of the community phylogenetic distance 
## matrix weighed by the abundances of the 
## species present in the community.
	weighted.mean(apply(my.com.dist, 1, min, na.rm = T), x[x > 0])
	
}


apply(my.sample, MARGIN = 1, Fmntd.a.sub)

}

Fmpd <-
function(dist.mat, my.sample){

Fmpd.sub = function(my.sub.sample){ 
		
## Get the names of species present in the 
## community
		com.names = names(my.sub.sample[my.sub.sample > 0])
		
## Calculate mpd by extracting the lower triangle 
## of a phylogenetic distance matrix comprised of 
## only the species in our community.
		mean(as.dist(dist.mat[com.names, com.names]))

	}


apply(my.sample, MARGIN = 1, Fmpd.sub)

}

Fmpd.a <-
function(dist.mat, my.sample){

Fmpd.a.sub = function(x){ 

## Get the names of the species in the community.
	com.names = names(x[x > 0])
		
## Make a matrix with one row containing 
## abundances and names of all present species
	com = t(as.matrix(x[x > 0]))
		
## Make phylogenetic distance matrix for taxa in ## community.
	com.dist = dist.mat[com.names, com.names]
		
## Calculate the product of the abundances of all 
## species in the community.
	abundance.products = t(as.matrix(com[1, com[1, ] > 0, drop = F]))%*% as.matrix(com[1, com[1, ] > 0, drop = F])
		
## Calculate a mean of the community phylogenetic 
## distance matrix weighted by the products of 
## all pairwise abundances.
	weighted.mean(com.dist, abundance.products)
	
	}


apply(my.sample, MARGIN = 1, Fmpd.a.sub)

}

Fsntd <-
function(dist.mat, my.sample){

Fsntd.sub = function(x){ 
		
## Get the names of the species present in a 
## community.
	com.names = names(x[x > 0])
		
## Make the community phylogenetic distance 
## matrix by extracting those rows and columns 
## that have species present in our community.
	my.com.dist = dist.mat[com.names, com.names]
		
## Set all diagonal values to NA so that the 
## zeros for conspecific comparisons do not 
## interfere with our calculation of nearest 
## neighbors.
	diag(my.com.dist) = NA
		
## Use apply() to calculate the minimum value in 
## each row of the community phylogenetic 
## distance matrix and take a mean of those 
## values.
	sd(apply(my.com.dist, MARGIN = 1, min,na.rm=T), na.rm=T)
	
}

apply(my.sample, MARGIN = 1, Fsntd.sub)

}

com.trait.moments <-
function(my.sample, traits){

mean.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
sd.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
skew.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
kurt.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))

skewness = function(x){ 
			m3<-sum((x-mean(x))^3)/length(x)
			s3<-sqrt(var(x))^3
			m3/s3
			}

kurtosis = function(x){
			m4 <- sum((x-mean(x))^4)/length(x)
			s4 <- var(x)^2
			m4/s4 -3


}



	for(i in 1:ncol(traits)){

		kurt.funk = function(x){ 
			kurtosis(traits[names(x[x > 0]), i])
		}

		skew.funk = function(x){ 
			skewness(traits[names(x[x >0]), i])	
		}

		mean.funk = function(x){ 
			mean(traits[names(x[x >0]), i], na.rm = T)
		}

		sd.funk = function(x){ 
			sd(traits[names(x[x >0]), i], na.rm = T)
		}


mean.output[,i] = apply(my.sample, MARGIN = 1, mean.funk)	
sd.output[,i] = apply(my.sample, MARGIN = 1, sd.funk)
skew.output[,i] = apply(my.sample, MARGIN = 1, skew.funk)
kurt.output[,i] = apply(my.sample, MARGIN = 1, kurt.funk)

	}

output = cbind(mean.output, sd.output, skew.output, kurt.output)

colnames(output)=paste(c(rep("mean",ncol(traits)), rep("sd",ncol(traits)), rep("skew",ncol(traits)),rep( "kurtosis",ncol(traits)) ), rep(names(traits), 4), sep=".")
rownames(output) = rownames(my.sample)
output



}

com.trait.weighted <-
function(my.sample, traits){

mean.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))
sd.output = matrix(NA, nrow=nrow(my.sample), ncol=ncol(traits))


	for(i in 1:ncol(traits)){
		weight.mean.funk = function(x){
			weighted.mean(traits[names(x[x>0]), i] , x[x > 0])
		}

		weight.sd.funk = function(x){
			wt.sd(traits[names(x[x>0]), i] , x[x>0])	
		}



		mean.output[,i] = apply(my.sample, MARGIN = 1, weight.mean.funk)	
		sd.output[,i] = apply(my.sample, MARGIN = 1, weight.sd.funk)
	}


output = cbind(mean.output, sd.output)
colnames(output)=paste(c(rep("mean",ncol(traits)), rep("sd",ncol(traits)) ), rep(names(traits), 2), sep=".")
rownames(output) = rownames(my.sample)
output

}

new.comdist = function(my.sample, my.dist.mat){

	get.presents = function(x){
		names(x[x>0])
		}

	list.of.names = apply(my.sample, 1, get.presents)
	
	Dpw.apply.function = function(x){
		tmp.function = function(z){
			mean(my.dist.mat[x,z])
			}
		
		lapply(list.of.names,FUN=tmp.function)
		
		}

dpw.output = lapply(list.of.names,Dpw.apply.function)

outt = do.call(rbind,lapply(dpw.output,unlist))


outt

}
new.comdist.prime = function(my.sample, my.dist.mat){

my.ra.sample = my.sample/rowSums(my.sample)

	get.presents = function(x){
		names(x[x>0])
		}

	get.weights = function(x){
		x[x>0]
		}

	list.of.weights = apply(my.ra.sample, 1, get.weights)

	
	Dpw.apply.function = function(x){
		tmp.function = function(z){
			dpw.output = sum(my.dist.mat[names(x),names(z)] * outer(x, z))

			
			}
		
		lapply(list.of.weights, FUN=tmp.function)
		
		}

dpw.output = lapply(list.of.weights, Dpw.apply.function)

outt = do.call(rbind,lapply(dpw.output,unlist))


outt
}




new.comdistnn = function(my.sample, my.dist.mat){

	get.presents = function(x){
		names(x[x>0])
		}

	list.of.names = apply(my.sample, 1, get.presents)
	
	Dnn.apply.function = function(x){
		tmp.function = function(z){
			
			mean(c(apply(my.dist.mat[x,z], MARGIN=1,min, na.rm=T), apply(my.dist.mat[x,z], MARGIN=2, min, na.rm=T)))
			
			}
		
		lapply(list.of.names, FUN=tmp.function)
		
		}

dnn.output = lapply(list.of.names,Dnn.apply.function)


outt = do.call(rbind,lapply(dnn.output,unlist))

outt

}
new.comdistnn.prime = function(my.sample, my.dist.mat){

my.ra.sample = my.sample/rowSums(my.sample)

	get.presents = function(x){
		names(x[x>0])
		}

	get.weights = function(x){
		x[x>0]
		}


	list.of.weights = apply(my.ra.sample, 1, get.weights)

	Dnn.apply.function = function(x){
		tmp.function = function(z){

			weighted.mean(c(apply(as.matrix(my.dist.mat[names(x),names(z)]), MARGIN=1,min, na.rm=T), apply(as.matrix(my.dist.mat[names(x),names(z)]), MARGIN=2, min, na.rm=T)),
			c(x,z)
			)		
			
			}
		
		lapply(list.of.weights, FUN=tmp.function)
		
		}

dnn.output = lapply(list.of.weights, Dnn.apply.function)

outt = do.call(rbind,lapply(dnn.output,unlist))


outt

}
rank.abund.pd = function(my.sample,my.phylo){

rank.com.fun = function(my.sub.sample){
my.species = sort(my.sub.sample[my.sub.sample > 0],decreasing=T)
pd.output = matrix(NA, nrow=1, ncol=length(my.species)-1)
	if(length(my.species)>=2){
	
		for(i in 2:length(my.species)){
	
			pd.output[,i-1] = pd(t(as.matrix(my.species[1:i])),my.phylo)$PD

		}

	}

	else{
		pd.output = c(NA)
	}
pd.output
}

outt= apply(my.sample, 1, rank.com.fun)
outt
}
sntd.a <-
function(dist.mat, my.sample){

sntd.a.sub = function(x){ 
		
## Get the names of species present in the 
## community.
	com.names = names(x[x > 0])
		
## Make the community phylogenetic distance 
## matrix using the names of the species present ## in the community.
	my.com.dist = dist.mat[com.names, com.names]
		
## Place NA values in the diagonals
	diag(my.com.dist) = NA
		
## Calculate a mean of the minimum values in each 
## row of the community phylogenetic distance 
## matrix weighed by the abundances of the 
## species present in the community.
	wt.sd(apply(my.com.dist, 1, min, na.rm = T), x[x > 0])
	
}

apply(my.sample, MARGIN = 1, sntd.a.sub)

}

trait.range <-
function(my.sample, traits){

	range.sub = function(x){

	com.names = names(x[x > 0])
	
	apply(traits[com.names, ], MARGIN = 2, max,na.rm=T) - apply(traits[com.names, ], MARGIN = 2, min,na.rm=T)
	}

apply(my.sample, MARGIN = 1, range.sub)

}

weighted.faith <-
function(my.phylo, my.sample){

weighted.faith.function = function(my.sub.sample){
		
			## extract the names of species in a community with an abundance greater than zero and use that information to make a pruned phylogeny for that community.
		tmp.tree = treedata(my.phylo, my.sub.sample[my.sub.sample > 0],warnings=F)$phy

	
			## Create empty branches matrix
		branches = matrix(NA, nrow = nrow(tmp.tree$edge), ncol = 4)

			## Fill first two columns of the matrix with node numbers defining each edge.
		branches[,1:2] = tmp.tree$edge

			## Fill the third column with the length of each branch
		branches[,3] = tmp.tree$edge.length
	
				get.leaves<-function(x){
					leaves.node<-tips(tmp.tree,x[2])
				}
		
				## Apply the get.leaves() function to each row in the branches object. This will retrieve species names subtended by each branch (i.e. ## row) in the branches matrix
		leaves = apply(branches, MARGIN = 1, get.leaves)
	
					## Now quickly loop through each set of leaves to ## calculate the mean abundance (Ai) for those species.
				for(i in 1:length(leaves)){
					branches[i,4] = mean(my.sub.sample[leaves[[i]]], na.rm = T) 

				}

			## Lastly calculated the Weighted Faith’s Index
		nrow(tmp.tree$edge) * ((sum(branches[,3] * branches[,4])) / sum(branches[,4]))


	}
outt = apply(my.sample, MARGIN = 1, weighted.faith.function)
outt

}

}, envir = lefse05)
