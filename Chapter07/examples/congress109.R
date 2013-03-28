countdata = read.csv("congress109.csv", header=TRUE, row.names=1)
memberdata = read.csv("congress109members.csv", header=TRUE, row.names=1)
countdata = as.matrix(countdata)

doublecenter = function(X)
{
	Z = X
	R = nrow(X)
	C = ncol(X)
	rowmu = rowMeans(X, na.rm=TRUE)
	colmu = colMeans(X, na.rm=TRUE)
	mu = mean(colmu, na.rm=TRUE)
	for(i in 1:C)
	{
		Z[,i] = (Z[,i] - colmu[i])
	}
	for(i in 1:R)
	{
		Z[i,] = (Z[i,] - rowmu[i])
	}
	Z = Z + mu
	return(Z)
}

standardize = function(X)
{
	Z = X
	C = ncol(X)
	colmu = colMeans(X, na.rm=TRUE)
	for(i in 1:C)
	{
		Z[,i] = Z[,i]/sd(Z[,i])
	}
	return(Z)
}

Z = doublecenter(countdata)
Z = standardize(Z)
pc1 = svd(Z)

comps = pc1$v
scores = pc1$u
plot(comps[,1:2])

scores = Z %*% comps

plot(scores[,1:2])
points(scores[which(memberdata$chamber=="H"),1:2], col='red', pch=19)
points(scores[which(memberdata$chamber=="S"),1:2], col='blue', pch=19)


pc2 = prcomp(countdata)
biplot(pc2)

