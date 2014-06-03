#############################################
#LECTURE 1
#############################################
#REJECTION ALGORITHM
#Beta Binomial

y = 9
n = 20

#Prior
theta <- runif(10000) 

#Sampling distribution
Y.sim <- rbinom(10000, 20, theta)

#Display the joint distribution
plot(Y.sim, theta, pch = 19, cex=.7, col = "grey4")
points(Y.sim[Y.sim == y], theta[Y.sim == y], cex=.7, pch = 19, col = "green4")

#Posterior distribution
theta.post <- theta[ Y.sim == y ]
summary(theta.post)

## Acceptance rate
mean(Y.sim == y)

## 95%CI
quantile(theta.post, prob=c(.025, .975))

par(mfrow = c(1,2))
u <- seq(0, 1, length= 200)

# Compare with theory
plot(u, dbeta(u, 10, 12), type="n", ylab="Frequency",
main="Posterior distribution")
hist(theta.post, col = "blue" , prob = T, xlab = "theta",add=T)
lines(u, dbeta(u, 10, 12), col = 2, lwd = 3)

# Posterior predictive distribution
y.post <- rbinom(length(theta.post), 20, theta.post )
hist(y.post, col = 5 , xlab = "y.rep", prob = T, main ="posterior predictive distribution")



#####################################################################
# Gaussian model
# Posterior distribution of (mu, sigma2)

n = 20
y = rnorm(n, 0, 1)

# Marginal distribution | y 
sigma.2 = (n-1)*var(y)/rchisq(10000, n-1)
# Conditional distribution | y
m = rnorm(10000, mean(y), sd = sqrt(sigma.2/n))

#2d sampling
plot(m, sigma.2, pch = 19, col = "orange", cex = .6)

#histograms
x11(); par(mfrow = c(1,2))
hist(m, prob = T, col = "green")
hist(sigma.2, prob = T, col = "orange")

###############
#Model checking
#Data from the cauchy model

skewness <- function(x) {
m3 <- mean((x-mean(x))^3)
skew <- m3/(sd(x)^3)
skew
} 

# Generate artificial data
y = rcauchy(n)

# CPosterior sampling
sigma.2 = (n-1)*var(y)/rchisq(10000, n-1)
mu = rnorm(10000, mean(y), sd = sqrt(sigma.2/n))
plot(mu, sigma.2, pch = 19, col = "orange", cex = .6)

par(mfrow = c(1,2))
hist(mu, prob = T, col = "green")
hist(sigma.2, prob = T, col = "orange")

## post predictive distribution of the skewness

post.pred = NULL
for (i in 1:1000)
{
ind = sample(10000, 1)
post.pred[i] = skewness(rnorm(n, m[ind], sqrt(sigma.2[ind])))
}

# show a histogram
hist(post.pred)
summary(post.pred)

# compare to the data 
skewness(y)
#[1] 2.365443



#####################################
# Application to the iris data set
# 

data(iris)
y = iris$Sepal.Length
n = length(y)

sigma.2 = (n-1)*var(y)/rchisq(10000, n-1)
mu = rnorm(10000, mean(y), sd = sqrt(sigma.2/n))
plot(mu, sigma.2, pch = 19, col = "orange", cex = .6)

par(mfrow = c(1,2))
hist(mu, prob = T, col = "green")
hist(sigma.2, prob = T, col = "orange")




#############################################
#LECTURE 2 MCMC
#############################################

#####################################################################
#Metropolis-Hastings beta binomial
#kernel = uniform = prior

n = 20 
y = 9

# initial value
theta.0 <- 0.9

# store posterior samples
theta.post <- theta.0

# monitor the log likelihood
log.likelihood <- NULL

for (i in 1:10000)
{
# proposal  
theta.1 <- runif(1)

# MH ratio 
ratio<-(theta.1/theta.0)^y *((1 - theta.1)/(1 - theta.0))^(n-y) 
if (runif(1) < ratio) {theta.0 <- theta.1} 

theta.post <- c(theta.post,theta.0)

log.likelihood <- c(log.likelihood, y*log(theta.0) + (n-y)*log(1 - theta.0) )
}

# Check that the algorithm is correct

par(mfrow = c(1,3))

u <- seq(0, 1, length= 200)
plot(u, dbeta(u, 10, 12), type="n", ylab="Frequency", main="Posterior distribution (All sweeps)")

hist(theta.post, col = "blue" , prob = T, add = T)
lines(u, dbeta(u, 10, 12), col = 2, lwd = 3)

# shows the history
plot(log.likelihood[1:200], cex = .5, pch = 19, type = "l")

plot(u, dbeta(u, 10, 12), type="n", ylab="Frequency", main="Posterior distribution (100 sweeps)")
hist(theta.post[1:100], col = "blue" , prob = T, add = T)
lines(u, dbeta(u, 10, 12), col = 2, lwd = 3)


#######################################################################
#Metropolis-Hastings
#kernel = non.uniform

nit = 10000
theta.0 <- 0.1
theta.post <- theta.0
log.likelihood <- NULL

for (i in 1:nit)
{
# proposal  
b = 100
theta.1 <- rbeta(1, b*theta.0/(1-theta.0), b)

# MH ratio 
ratio1 <- (theta.1/theta.0)^(y)* 
( (1 - theta.1)/(1 - theta.0) )^(n-y) 
ratio2 <-  dbeta(theta.0, b*theta.1/(1-theta.1), b)/
dbeta(theta.1, b*theta.0/(1-theta.0), b) 
ratio <- ratio1*ratio2

if (runif(1) < ratio) {theta.0 <- theta.1} 

theta.post <- c(theta.post,theta.0)
log.likelihood <- c(log.likelihood, y*log(theta.0) + (n-y)*log(1 - theta.0) )
}

# show results

u <- seq(0, 1, length= 200)

par(mfrow=c(1,3))

plot(u, dbeta(u, 10, 12), type="n", ylab="Frequency", main="Posterior distribution (burnin = 200)", xlab = "x")
hist(theta.post[200:5000], col = "blue" , prob = T, add = T)
lines(u, dbeta(u, 10, 12), col = 2, lwd = 3)

plot(u, dbeta(u, 10, 12), type="n", ylab="Frequency", main="150 sweeps", xlab = "x")
hist(theta.post[1:150], col = "yellow" , prob = T, add = T)
lines(u, dbeta(u, 10, 12), col = 2, lwd = 3)


plot(log.likelihood[1:100], cex = .5, pch = 19, type = "l")



#######################################################################
## Gibbs sampling 
## Gaussian example (vary rho) 

par(mfrow = c(1, 2))

rho <- .99

## Exact simulation
theta1 <- rnorm(200)
theta2 <- rnorm(200, rho*theta1, sd = sqrt(1 - rho^2))
plot(theta1, theta2)
cor(theta1, theta2)


## Gibbs sampler is mixing slowly

nit = 2000
theta1 <- -3; theta2 <- 3

theta1.post <- theta1
theta2.post <- theta2

for (i in 1:nit){

theta1 <- rnorm(1, rho*theta2, sd = sqrt(1 - rho^2))
theta2 <- rnorm(1, rho*theta1, sd = sqrt(1 - rho^2))

theta1.post <- c(theta1.post, theta1)
theta2.post <- c(theta2.post, theta2)
}

# show results

plot(theta1.post[-1], theta2.post[-1])
cor(theta1.post, theta2.post)

hist(theta2.post)


#######################################################################
## Gibbs sampler (m, sigma2)
##  

y = rnorm(200)
post.mu = 10
post.sigma.2 = 0.3

for (i in 2:2000){
s.2 = sum( (y - post.mu[i-1])^2 )
post.sigma.2 = c(post.sigma.2, s.2/rchisq(1, n))
post.mu = c(post.mu,rnorm(1, mean(y), sd = sqrt(post.sigma.2[i]/n)))}


# Show results

par(mfrow = c(1,2))
hist(post.mu[-1], prob = T, col = "green", main = "mean")
hist(post.sigma.2[-(1:100)], prob = T, col = "orange", main = "variance")



#############################################
#LECTURE 3 Gaussian mixture 
#############################################


z = c(rep(0,60), rep(1, 140))
y = NULL

# within-group means and sd
m.t = c(-1, 2)
s.t = c(.7, 1.6) 

for (i in 1:200){
y[i] = rnorm(1, m.t[z[i]+1], sd = s.t[z[i]+1])}


## Show it

x = seq(-5, 7, length = 500)

par(mfrow = c(1, 2))
hist(y, prob = T, col = "grey")
hist(y, prob = T, col = "grey")
points(x, .3*dnorm(x, -1, .7) , lwd=4, type ="l", col = "yellow")
points(x, .7*dnorm(x, 2, 1.6) , lwd=4, type ="l", col = "blue2")
points(x, .3*dnorm(x, -1, .7) + .7*dnorm(x, 2, 1.6) , lwd=4, type ="l", col = "green4")


## MCMC for mixture 

mcmc.mix = function(y,niter=1000, m.o= c(0,1), sigma2.o = c(1,1)){
n = length(y)
K = length(m.o)

m = m.o ;s2 = sigma2.o; nk = NULL 
z = NULL
p = NULL 

p.mcmc = NULL
m.mcmc = NULL
s2.mcmc = NULL 
z.mcmc = NULL

for (nit in 1:niter){

for (i in 1:n){
p = exp(-(y[i] - m)^2/2/s2 )/sqrt(s2)
z[i] = sample(1:K, 1, prob = p )}

z.mcmc = rbind(z.mcmc, z)

for (k in 1:K){
nk[k] = sum(z==k) 
m[k] = rnorm( 1, mean(y[z==k]), sd = sqrt(s2[k]/nk[k]) )
s2[k] = sum((y[z==k] - m[k])^2)/rchisq(1, nk[k])
}

p.mcmc = rbind(p.mcmc, nk/n)
m.mcmc = rbind(m.mcmc, m)
s2.mcmc = rbind(s2.mcmc, s2)
}
return(list(z= z.mcmc, p = p.mcmc, m = m.mcmc, sigma2 = s2.mcmc))
}

### modele K = 2

obj = mcmc.mix(y, niter = 300, m.o = c(2, 4), sigma2.o = c(1,1))

mat = rbind(apply(obj$z[-(1:100),], MARGIN = 2, FUN = function(a) mean(a == 1) ), apply(obj$z[-(1:100),], MARGIN = 2, FUN = function(a) mean(a == 2)),apply(obj$z[-(1:100),], MARGIN = 2, FUN = function(a) mean(a == 3) ))

barplot(mat, col = c("yellow","blue2","red"))


### modele K = 3
x11()

obj = mcmc.mix(y, niter = 200, m.o = c(-1, 1, 2), sigma2.o = c(1,1,1))

barplot(rbind(apply(obj$z[-(1:100),], MARGIN = 2, FUN = function(a) mean(a == 1) ), apply(obj$z[-(1:100),], MARGIN = 2, FUN = function(a) mean(a == 2)),apply(obj$z[-(1:100),], MARGIN = 2, FUN = function(a) mean(a == 3) )), col = c("yellow","blue2","red"))


### Posterior checking


obj = mcmc.mix(y, niter = 300, m.o = c(-1, 1, 2), sigma2.o = c(1,1,1))


post.stat = NULL
for (i in 1:100){
zz = obj$z[i+100,]
mm = obj$m[i+100,]
ss2 = obj$sigma2[i+100,] 
y.rep = c(rnorm(sum(zz==1), mm[1], sqrt(ss2[1])), 
       rnorm(sum(zz==2), mm[2], sqrt(ss2[2])), 
       rnorm(sum(zz==3), mm[3], sqrt(ss2[3])))
post.stat = c(post.stat, skewness(y.rep))
 }


hist(post.stat)
skewness(y)




