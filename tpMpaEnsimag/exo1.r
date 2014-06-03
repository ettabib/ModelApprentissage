Y <- c(9,5,7,15,9,0,22,1,9,11,9,13,12,11,2,1,26,7,0)
n=length(Y)
theta.0 <- 0.1
theta.post <- theta.0
Lambda=1/8
Nbite=10^4
m=mean(Y)
for( i in 1:Nbite)
{
  ## proposal
  theta.1 <- rexp(1,Lambda)
  ## calcul ration
  ratio1 <- exp(-(n-Lambda)*(theta.1-theta.0))
  ratio2 <- (theta.1/theta.0)^(n*m-1)
  ratio <-  ratio1 * ratio2
  if (runif(1)<ratio){
    theta.0 <- theta.1
  }
  theta.post <- c(theta.post,theta.0)
  }

## par(mfrow= c(1,2));
##
## par(new=TRUE)
u <- seq(0, 20, length= 200);
plot(u, dgamma(u,shape=n*mean(Y),scale=1/n),  ylab="Frequency",xlim=c(6,12));
hist(theta.post,xlim=c(6,12),prob=T,add=T);

## intervalle de credibilité à 95%
interval.cred <- c(sort(theta.post)[25],sort(theta.post)[975])

