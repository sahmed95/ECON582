# Problem 1: Exercise 9.26

data <- read.table("/Users/student/Desktop/Spring21/582/HW1/Nerlove1963.txt", header = TRUE)
log_C <- matrix(log(data$Cost), ncol=1)
log_Q <- matrix(log(data$output), ncol=1)
log_PL <- matrix(log(data$Plabor), ncol=1)
log_PK <- matrix(log(data$Pcapital), ncol=1)
log_PF <- matrix(log(data$Pfuel), ncol=1)
x <- as.matrix(cbind(matrix(1,nrow(log_C),1),log_Q, log_PL, log_PK, log_PF))
y<- log_C
n<- nrow(x)
k<- ncol(x)

# a) Unrestricted regression 

invx <- solve(t(x)%*%x)
b_ols <- solve(t(x)%*%x)%*%(t(x)%*%y)
e_ols <- rep((y-x%*%b_ols), times=k)
xe_ols <-x*e_ols 
V_ols <- (n/(n-k))*invx%*%(t(xe_ols)%*%xe_ols)%*%invx
se_ols <- sqrt(diag(V_ols))
print(b_ols)
print(se_ols)

#-------------------------------------------------------------

# c) Constrained regression 

R <- c(0,0,1,1,1)
c <- 1
iR <- invx%*%R%*%solve(t(R)%*%invx%*%R)
b_cls <- b_ols -iR%*%t(R)%*%b_ols +iR%*%c
e_cls <- rep((y-x%*%b_cls),times=k)
xe_cls <- x*e_cls 
V_tilde <- (n/(n-k+1))*invx%*%(t(xe_cls)%*%xe_cls)%*%invx
V_cls <- V_tilde-iR%*%t(R)%*%V_tilde-V_tilde%*%t(iR%*%t(R))+iR%*%t(R)%*%V_tilde%*%t(iR%*%t(R))
se_cls <- sqrt(diag(V_cls))
print(b_cls)
print(se_cls)

#-------------------------------------------------------------

# d) Efficient minimum distance 
iV <- solve(t(R)%*%V_ols%*%R)
V<-V_ols%*%R%*%iV
b_emd <- b_ols-V%*%t(R)%*%b_ols+V%*%c
e_emd<- rep((y-x%*%b_emd), times=k)
xe_emd <-x*e_emd
V2 <- (n/(n-k+1))*invx%*%(t(xe_emd)%*%xe_emd)%*%invx
V_emd <-V2 - V2%*%R%*%solve(t(R)%*%V2%*%R)%*%t(R)%*%V2
se_emd<- sqrt(diag(V_emd))
print(b_emd)
print(se_emd)

#-------------------------------------------------------------

# e) Wald statistic
q <- length(c)
V_r <- solve(t(R)%*%V_ols%*%R)
W <- t(t(R)%*%b_ols-c)%*%V_r%*%t(t(R)%*%b_ols-c)
alpha <- 0.05 
C <- qchisq(1-alpha, q) 
print(W)
print(C)

if (W>C){
	print("Reject H0")
} else {
	print("Accept H0")
}

#-------------------------------------------------------------

# f) Minimun distance statistic 
q <- length(c)
J <- t((b_ols-b_emd))%*%solve(V_ols)%*%(b_ols-b_emd)
alpha <- 0.05 
C <- qchisq(1-alpha, q) 
print(J)
print(C)
if (J>C){
	print("Reject H0")
} else {
	print("Accept H0")
}


#-------------------------------------------------------------

# Exercise 10.28
set.seed(1)

# Creating a function to calculate Beta hat

estimate.beta <- function(z,y){
	n<- nrow(z)
	k<- ncol(z)


	invz <- solve(t(z)%*%z)
	b_ols <- solve(t(z)%*%z)%*%(t(z)%*%y)
	e_ols <- rep((y-z%*%b_ols), times=k)
	ze_ols <-z*e_ols 
	V_ols <- (n/(n-k))*invz%*%(t(ze_ols)%*%ze_ols)%*%invz
	se_ols <- sqrt(diag(V_ols))
	return(b_ols)
	return(se_ols)
}

# Asymptotic

beta.hat <- estimate.beta(z=x, y=log_C)
print(beta.hat)
print(se_ols)


# Jackknife

beta.hat.jack <- matrix(data=0, nrow = k, ncol = n)

for (i in 1:n) {
	df.loo.i <- x[-i,]
	p<- log_C[-i,]
	beta.hat.jack[,i] <- estimate.beta(z=df.loo.i,y=p)
}

beta.bar.jack <- rowMeans(beta.hat.jack)
diff.jack <- (beta.hat.jack-beta.bar.jack)
var.jack <- ((n-1)/n)*(diff.jack%*%t(diff.jack))
se.jack <- sqrt(diag(var.jack))
se.jack

# Bootstrap 

B<- 10000

beta.hat.boot <- matrix(0, nrow=k, ncol =B)

for (b in 1:B){
	# Construct a b-th bootstrap sample 
	idx.b <- sample(n, replace = TRUE)
	df.b <- x[idx.b,]
	y.b <- log_C[idx.b, ]
	beta.hat.boot[,b] <- estimate.beta(z=df.b, y=y.b)
}

beta.bar.boot <- rowMeans(beta.hat.boot)
diff.boot <- (beta.hat.boot-beta.bar.boot)
var.boot <- (1/(B-1))*(diff.boot %*% t(diff.boot))
se.boot <- sqrt(diag(var.boot))
se.boot

#-------------------------------------------------------------

# b) 

theta.hat <- b_ols[3]+b_ols[4]+b_ols[5]
theta.hat

# Asymptotic (Using the Delta Method)

var.theta <- t(R)%*%V_ols%*%R
se.theta <- sqrt(var.theta)
print(se.theta)

# Jackknife 

theta.hat.jack <-rep(0,n)

for (i in 1:n) {
	theta.hat.jack[i] <- beta.hat.jack[3,i]+beta.hat.jack[4,i]+beta.hat.jack[5,i]
}

theta.bar.jack <- mean(theta.hat.jack)
var.jack <- (n-1)*mean((theta.hat.jack-theta.bar.jack)^2)
se.jack <- sqrt(var.jack)
se.jack



# Bootstrap 

B<- 10000
theta.hat.boot <- rep(0,B)
for (b in 1:B){
	theta.hat.boot[b] <- beta.hat.boot[3,b]+beta.hat.boot[4,b]+beta.hat.boot[5,b]
}

theta.bar.boot <- mean(theta.hat.boot)
var.boot <- (B/(B-1))*mean((theta.hat.boot-theta.bar.boot)^2)
se.boot <- sqrt(var.boot)
se.boot

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



