# Problem 4: Exercise 4.26
library(tidyverse) # Packages for data manipulation 
library(readxl) # Read xlsx file 
library(sandwich) # Robust standard error 
library(lmtest)


DDK2011 <- read_excel("DDK2011.xlsx")
n.DDK2011 <- nrow(DDK2011)
 
DDK2011.sub <- DDK2011 %>% 
	mutate(testscore = scale(totalscore)) %>%
	# Removing observations with missing variables
	filter(testscore != ".", tracking != ".", agetest != ".", etpteacher != ".", percentile != ".", girl != '.') %>% 
	mutate_all(as.numeric)


OLS.out <- lm(testscore ~ tracking+agetest+girl+etpteacher+percentile, data = DDK2011.sub)
b.hat <-coef(OLS.out)
b.hat 

# Robust standard error 

V.HC1 <- vcovHC(OLS.out, type = "HC1")
se.HC1 <- sqrt(diag(V.HC1))
se.HC1

# Clustered standard error

V.cluster <- vcovCL(OLS.out, cluster =~ schoolid)
se.cluster<- sqrt(diag(V.cluster))
se.cluster

# Calculating the absolute difference between the two standard errors

diff <- abs(se.HC1-se.cluster)
diff

# Calculating the percentage difference in the two standard errors
percent<- (se.cluster-se.HC1)/(se.HC1)*100
percent
#---------------------------------------------------------------------------------------------

# Exercise 10.31 
DDK2011 <- read_excel("DDK2011.xlsx")
DDK2011 <- DDK2011 %>% 
	mutate(testscore = scale(totalscore)) %>%
	# Removing observations with missing variables
	filter(testscore != ".", tracking != ".", agetest != ".", etpteacher != ".", percentile != ".", girl != '.') %>% 
	mutate_all(as.numeric)
	
n.DDK2011 <- nrow(DDK2011)
DDK2011_group <- DDK2011 %>% group_nest(schoolid)

cluster <- unique(DDK2011_group$schoolid)
n.cluster <- length(cluster)

# Jackknife 

n.beta <- 6
beta.cluster.loo <- matrix(0, nrow =n.beta, ncol =n.cluster)

for (i in 1:n.cluster){
	df.cluster.loo.i <- DDK2011_group[-i,] %>% unnest(data)
	ols.out.i <- lm(testscore ~ tracking+agetest+girl+etpteacher+percentile, data = df.cluster.loo.i)
	beta.loo.i <- coef(ols.out.i)
	
	beta.cluster.loo[,i] <- beta.loo.i
}

beta.bar.cluster <- rowMeans(beta.cluster.loo)
diff.jack.cluster <- beta.cluster.loo-beta.bar.cluster
var.cluster.jack <- (n.cluster-1)/n.cluster*diff.jack.cluster%*%t(diff.jack.cluster)
se.cluster.jack <- sqrt(diag(var.cluster.jack))
se.cluster.jack

#--------------------------------------------------------------------------------------------
# Bootstrap 

set.seed(1)
DDK2011 <- read_excel("DDK2011.xlsx")
DDK2011 <- DDK2011 %>% 
	mutate(testscore = scale(totalscore)) %>%
	# Removing observations with missing variables
	filter(testscore != ".", tracking != ".", agetest != ".", etpteacher != ".", percentile != ".", girl != '.') %>% 
	mutate_all(as.numeric)
	
n.DDK2011 <- nrow(DDK2011)
DDK2011_group <- DDK2011 %>% group_nest(schoolid)

cluster <- unique(DDK2011_group$schoolid)
n.cluster <- length(cluster)
n.beta <- 6

B <- 10000
beta.hat.boot.cluster <- matrix(0, nrow = n.beta, ncol = B)

for (b in 1:B){
	idx.boot <- sample(1:n.cluster, size = n.cluster, replace = TRUE)
	df.boot.b <- DDK2011_group[idx.boot, ] %>% unnest(data)
	ols.out.b <- lm(testscore ~ tracking+agetest+girl+etpteacher+percentile, data = df.boot.b)
	beta.boot.b <- coef(ols.out.b)
	beta.hat.boot.cluster[, b] <- beta.boot.b
}

beta.bar.boot.cluster <- rowMeans(beta.hat.boot.cluster)
diff.boot.cluster <- (beta.hat.boot.cluster - beta.bar.boot.cluster)
var.boot.cluster <- (1/(B-1))*(diff.boot.cluster%*%t(diff.boot.cluster))
se.boot.cluster <- sqrt(diag(var.boot.cluster))
se.boot.cluster

# Constructing BC_a percentile intervals for each coefficient 

alpha <- 0.05 
alphas <- c(alpha/2, 1-alpha/2)

# evaluate z_alpha for each alpha values 
z.alphas <- qnorm(alphas)


for (i in 1:n.beta){
	p.star.i <- mean(beta.hat.boot.cluster[i,] <= b.hat[i])
	z0.star.i <- qnorm(p.star.i)
	beta.bar.loo <- mean(beta.cluster.loo[i,])
	diff.jack <- beta.bar.loo-beta.cluster.loo[i,]
	a.jack.num <- sum(diff.jack^3)
	a.jack.den <- 6*sum(diff.jack^2)^(3/2)
	a.jack <- a.jack.num/a.jack.den
	correction <- (z.alphas+z0.star.i)/(1-a.jack*(z.alphas+z0.star.i))
	x.alphas <-pnorm(z0.star.i+correction)
	q.star.x.alphas <- quantile(beta.hat.boot.cluster[i,], probs =x.alphas)
	CI_BCa <- q.star.x.alphas
	print(CI_BCa)
}