# Problem 2: Exercise 9.27

mrw <- read.table("/Users/student/Desktop/Spring21/582/HW1/MRW1992.txt", header=TRUE)
N <- matrix(mrw$N, ncol =1)
lndY <- matrix(log(mrw$Y85)-log(mrw$Y60),ncol=1)
lnY60 <- matrix(log(mrw$Y60), ncol=1)
lnI <- matrix(log(mrw$invest/100), ncol =1)
lnG <- matrix(log(mrw$pop_growth/100+0.05), ncol=1)
lnS <- matrix(log(mrw$school/100), ncol =1)
X <- as.matrix(cbind(lnY60, lnI, lnG, lnS,matrix(1,nrow(lndY),1)))
x <- X[N==1,]
y <- lndY[N==1]

# Creating a function for estimating beta 

estimate.beta <- function(z,y){
	n<- nrow(z)
	k<- ncol(z)
	invz <- solve(t(z)%*%z)
	b_ols <- solve(t(z)%*%z)%*%(t(z)%*%y)
	return(b_ols)
	
}


# Unrestricted regression 

beta.hat <- estimate.beta(x,y)

# Standard error 
n <-nrow(x)
k <- ncol(x)
invx <- solve(t(x)%*%x)
e_ols <- rep((y-x%*%beta.hat), times=k)
xe_ols <-x*e_ols 
V_ols <- (n/(n-k))*invx%*%(t(xe_ols)%*%xe_ols)%*%invx
se_ols <- sqrt(diag(V_ols))
print(beta.hat)
print(se_ols)


# Test (Wald Statistic)

R <- c(0,1,1,1,0)
c<- 0
q <- length(c)
V_r <- solve(t(R)%*%V_ols%*%R)
W <- t(t(R)%*%beta.hat-c)%*%V_r%*%t(t(R)%*%beta.hat-c)
alpha <- 0.05 
C <- qchisq(1-alpha, q) 
print(W)
print(C)

if (W>C){
	print("Reject H0")
} else {
	print("Accept H0")
}
#-----------------------------------------------------------------

# Exercise 10.29 

# Part a) 

# Asymptotic 
n <-nrow(x)
k <- ncol(x)
beta.hat <- estimate.beta(x,y)
invx <- solve(t(x)%*%x)
e_ols <- rep((y-x%*%beta.hat), times=k)
xe_ols <-x*e_ols 
V_ols <- (n/(n-k))*invx%*%(t(xe_ols)%*%xe_ols)%*%invx
se_ols <- sqrt(diag(V_ols))
print(beta.hat)
print(se_ols)

# Jackknife 

beta.hat.jack <- matrix(data=0, nrow = k, ncol = n)

for (i in 1:n) {
	df.loo.i <- x[-i,]
	p<- y[-i]
	beta.hat.jack[,i] <- estimate.beta(z=df.loo.i,y=p)
}

beta.bar.jack <- rowMeans(beta.hat.jack)
diff.jack <- (beta.hat.jack-beta.bar.jack)
var.jack <- ((n-1)/n)*(diff.jack%*%t(diff.jack))
se.beta.jack <- sqrt(diag(var.jack))
se.beta.jack

# Bootstrap 

B<- 10000
set.seed(1)
beta.hat.boot <- matrix(0, nrow=k, ncol =B)

for (b in 1:B){
	# Construct a b-th bootstrap sample 
	idx.b <- sample(n, replace = TRUE)
	df.b <- x[idx.b,]
	y.b <- y[idx.b]
	beta.hat.boot[,b] <- estimate.beta(z=df.b, y=y.b)
}

beta.bar.boot <- rowMeans(beta.hat.boot)
diff.boot <- (beta.hat.boot-beta.bar.boot)
var.boot <- (1/(B-1))*(diff.boot %*% t(diff.boot))
se.beta.boot <- sqrt(diag(var.boot))
se.beta.boot

#-----------------------------------------------------------------

# Part b) Estimating theta 

# Asymptotic 

theta.hat <- beta.hat[2]+beta.hat[3]+beta.hat[4]
theta.hat

# Standard error (Using Delta Method)
var.theta <- t(R)%*%V_ols%*%R
se.theta <- sqrt(var.theta)
print(se.theta)

# Jackknife 

theta.hat.jack <-rep(0,n)

for (i in 1:n) {
	theta.hat.jack[i] <- beta.hat.jack[2,i] + beta.hat.jack[3,i]+beta.hat.jack[4,i]
}

theta.bar.jack <- mean(theta.hat.jack)
var.theta.jack <- (n-1)*mean((theta.hat.jack-theta.bar.jack)^2)
se.theta.jack <- sqrt(var.theta.jack)
se.theta.jack

B<- 10000

theta.hat.boot <- rep(0,B)
for (b in 1:B){
	theta.hat.boot[b] <- beta.hat.boot[2,b]+beta.hat.boot[3,b]+beta.hat.boot[4,b]
	}

theta.bar.boot <- mean(theta.hat.boot)
var.theta.boot <- (B/(B-1))*mean((theta.hat.boot-theta.bar.boot)^2)
se.theta.boot <- sqrt(var.theta.boot)
se.theta.boot

#-------------------------------------------------------------

# c) Confidence intervals for theta 

# Percentile method 

alpha <- 0.05
q.star.alphas <- quantile(theta.hat.boot, probs = c(alpha/2, 1-alpha/2))

CI_percentile <- q.star.alphas 
CI_percentile 

# BC_a Percentile interval 

alphas <- c(alpha/2, 1-alpha/2)

# Evaluate z_alpha for each alpha values 

z.alphas <- qnorm(alphas)

# Evaluate z0.star 

p.star <- mean(theta.hat.boot <= theta.hat)
z0.star <- qnorm(p.star)

diff.jack.t <- theta.bar.jack - theta.hat.jack 
a.jack.num <- sum(diff.jack^3)
a.jack.den <- 6*sum(diff.jack^2)^(3/2)
a.jack <- a.jack.num/a.jack.den

correction <- (z.alphas+z0.star)/(1-a.jack*(z.alphas+z0.star))
x.alphas <-pnorm(z0.star+correction)
q.star.x.alphas <- quantile(theta.hat.boot, probs = x.alphas)

CI_BCa <- q.star.x.alphas
CI_BCa

