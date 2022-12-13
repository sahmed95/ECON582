# Problem 3: Exercise 7.28 

library(tidyverse) # Packages for data manipulation 
library(readxl) # Read xlsx file 
library(sandwich) # Robust standard error 
library(lmtest)

cps09mar <- read_excel("cps09mar.xlsx")

cps09mar.sub <- cps09mar %>%
	mutate(log_wage = log(earnings/(hours*week))) %>%
	mutate(experience = age-education-6) %>%
	mutate(experience2 = experience^2/100) %>%
	filter(female == 0, hisp ==1,race ==1)
	
OLS.out <- lm(log_wage ~ education+experience+experience2, data = cps09mar.sub)
coeftest(OLS.out, vcov = vcovHC(OLS.out, type="HC1"))
beta.hat <-coef(OLS.out)
beta.hat

#------------------------------------------------------------------------------------
# b) 

theta.hat <- 5*beta.hat[2]/(5*beta.hat[3]+beta.hat[4])
theta.hat
#------------------------------------------------------------------------------------
# c) 

V_ols <- vcovHC(OLS.out, type = "HC1")
R.hat <- c(0, 5/(5*beta.hat[3]+beta.hat[4]), (-25*beta.hat[2])/(5*beta.hat[3]+beta.hat[4])^2
, (-5*beta.hat[2])/(5*beta.hat[3]+beta.hat[4])^2)
V.theta <- t(R.hat)%*%V_ols%*%R.hat
se.theta.hat <-sqrt(V.theta)
se.theta.hat

#------------------------------------------------------------------------------------
# d) 
alpha <- 0.1
c<- qnorm(1-alpha/2)
u.bound <- theta.hat+c*se.theta.hat
l.bound <- theta.hat-c*se.theta.hat
text <- 'Confidence interval:[ ${l.bound} , ${u.bound}]'
cat(str_interp(text))

#------------------------------------------------------------------------------------
# e) 
educ <- 12
exp <- 20

coef <- c(1,educ, exp, exp^2/100)
reg <- t(coef)%*%beta.hat

# SE (by the Delta method)

V.reg<- t(coef)%*%V_ols%*%coef
se.reg <- sqrt(V.reg)

# Confidence interval 
alpha <- 0.05
c <- qnorm(1-alpha/2)
u.bound <- reg+c*se.reg
l.bound <- reg-c*se.reg
text <- 'Confidence interval:[ ${l.bound} , ${u.bound}]'
cat(str_interp(text))

#------------------------------------------------------------------------------------

# f) 
# Computing the regression function with the new values 

out.educ <- 16
out.exp <- 5
x <- c(1, out.educ, out.exp, (out.exp)^2/100)
f_reg <- t(x)%*%beta.hat

# Standard error (we need sigma.hat^2)

e<- residuals(OLS.out)
sigma_2 <- mean(e^2)
V.x <- sigma_2+t(x)%*%V_ols%*%x
se.x <- sqrt(V.x)

# 80% confidence interval for log(wage)

alpha <- 0.2
c<- qnorm(1-alpha/2)
u.bound <- f_reg+c*se.x
l.bound <- f_reg-c*se.x
text <- 'Confidence interval:[ ${l.bound} , ${u.bound}]'
cat(str_interp(text))

# 80% confidence interval for wage 

exp.u.bound<- exp(u.bound)
exp.l.bound<- exp(l.bound)
text <- 'Confidence interval:[ ${exp.l.bound} , ${exp.u.bound}]'
cat(str_interp(text))


#------------------------------------------------------------------------------------

# Exercise 10.30

cps09mar.sub1 <- cps09mar.sub %>%
	filter(region ==2, marital == 7)
	

# Creating a function to estimate theta 

estimate.theta <-function(data){
	OLS1.out <- lm(log_wage ~ education+experience+experience2, data = data)
	b_ols <- coef(OLS1.out)
	theta.hat <- 5*b_ols[2]/(5*b_ols[3]+b_ols[4])
	return(theta.hat)
}

# Asymptotic 

OLS2.out <- lm(log_wage ~ education+experience+experience2, data = cps09mar.sub1)
b_as <- coef(OLS2.out)
theta.as <- 5*b_as[2]/(5*b_as[3]+b_as[4])
theta.as
V_as <- vcovHC(OLS2.out, type = "HC1")
R.hat.as <- c(0, 5/(5*b_as[3]+b_as[4]), (-25*b_as[2])/(5*b_as[3]+b_as[4])^2
, (-5*b_as[2])/(5*b_as[3]+b_as[4])^2)
V.theta.as <- t(R.hat.as)%*%V_as%*%R.hat.as
se.theta.as <-sqrt(V.theta.as)
se.theta.as

# Jackknife 

n <- nrow(cps09mar.sub1)
theta.hat.jack <-rep(0,n)

for (i in 1:n) {
	df.loo.i <- cps09mar.sub1[-i,]
	theta.hat.jack[i] <- estimate.theta(data = df.loo.i)
}

theta.bar.jack <- mean(theta.hat.jack)
var.jack <- (n-1)*mean((theta.hat.jack-theta.bar.jack)^2)
se.jack <- sqrt(var.jack)
se.jack


# Bootstrap 

set.seed(1000)

B<- 10000
theta.hat.boot <- rep(0,B)
for (b in 1:B){
	# Construct a b-th bootstrap sample 
	idx.b <- sample(n, replace =TRUE)
	df.b <- cps09mar.sub1[idx.b,]
	theta.hat.boot[b]<-estimate.theta(data=df.b)	
}

theta.bar.boot <- mean(theta.hat.boot)
var.boot <- (B/(B-1))*mean((theta.hat.boot-theta.bar.boot)^2)
se.boot <- sqrt(var.boot)
se.boot


#------------------------------------------------------------------------------------

# BC Percentile Interval 

alpha <- 0.05
alphas <- c(alpha/2, 1-alpha/2)

# Evaluate z_alpha for each alpha values 
z.alphas<-qnorm(alphas)

# Evaluate z0.star

p.star <- mean(theta.hat.boot <= theta.hat)
z0.star <- qnorm(p.star)

# Calculate x(alpha) for each alpha values 
x.alphas <- pnorm(z.alphas+2*z0.star)

q.star.x.alphas <- quantile(theta.hat.boot, probs = x.alphas)
CI_BC<- q.star.x.alphas
CI_BC

