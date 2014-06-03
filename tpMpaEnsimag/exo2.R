
# def des fonctions de calcule
# retourne un estimateur de alpha de la loi beta(alpha,beta) 
estim.alpha <- function(x)
{
  return (x[1]*((x[1]*(1-x[1])/(x[2]))-1))
}


# retourne un estimateur de beta de la loi beta(alpha,beta)
estim.beta <- function(x)
{
  return ((1-x[1])*((x[1]*(1-x[1])/(x[2]))-1))
}


# lecture du fichier de données
  data <- read.table("r.asc", sep=" ")

# def des paramettres du pb
Nbite=1
theta.0 <- matrix(nrow=1,ncol=2)
theta.1 <- matrix(nrow=1,ncol=2)
theta.post <- matrix(nrow=Nbite+1,ncol=2)
Lambda <- matrix(nrow=1,ncol=2)
m1=0
m2=0
  
  
# initialisation des paramettres
yi<-data[,1]
ni<-data[,2]
n<-length(data[,1])
theta.exper<-c(1,1)
Niter=10^4
thetai.post <- matrix(nrow=Niter+1,ncol=n)
  
  
# theta apartir des données
  for ( i in (1:n)){
    theta.exper[i] <-  yi[i]/ni[i]
  }  

# Generation de alpha et de beta
    alpha.estim <-estim.alpha(c(mean(theta.exper),var(theta.exper)))  
    beta.estim<- estim.beta(c(mean(theta.exper),var(theta.exper)))  

  # Valeur initiale arbitraire
  for ( i in (1:n)){
   thetai.post[1,i] <- 0.01 
  }  
  
  
  for ( j in (1:Niter)) {

  # Generation de alpha et de beta
  
    # initialisation par une valeur arbitraire
    theta.0 <-c(1,1)    
    Lambda <- c(1/alpha.estim,1/beta.estim)
    theta.post[1, ] <- theta.0
  
    for( i in 1:Nbite)
    {
     
      #proposal
      theta.1 <-   c(rexp(1,Lambda[1]),rexp(1,Lambda[2])) # c(runif(1,0,100),runif(1,0,100)) #
      ## calcul ration
      ratio1 <- ((theta.1[1]+theta.1[2])/(theta.0[1]+theta.0[2]))^(-5/2)
      ratio2 <- exp(-Lambda[1]*(theta.0[1]-theta.1[1])-Lambda[2]*(theta.0[2]-theta.1[2]))
      ratio3 <-  (prod(thetai.post[j,]))^(theta.1[1]-theta.0[1])*(prod(1-thetai.post[j,]))^(theta.1[2]-theta.0[2])
      ratio <- ratio1*ratio2*ratio3 
      if (runif(1)<ratio){
        theta.0 <- theta.1
      }
      theta.post[i+1, ] <- theta.0
    }
    
    # estimation de aplha et beta
    m<-mean(theta.post[,1])/(mean(theta.post[,1])+mean(theta.post[,2]))
    m1<-c(m1,theta.post[2,1])
    m2<-c(m2,theta.post[2,2])
    
    alpha1<-theta.post[2,1]
    beta1<- theta.post[2,2]
    
    
    #alpha1<-c(alpha1,mean(theta.post[,1]))
    #beta1<-c(alpha2,mean(theta.post[,2]))
  
    # Echantillonage de Gibs pour les theta_i
    thetai.post[j+1,]=rbeta(n,shape1=yi+alpha1,shape2=ni- yi + beta1 )
  }
    
  
  
  moyen.theta.post <- matrix(nrow=1,ncol=n)
  
  p = 300
  for (k in (1:n) )
  {
  moyen.theta.post[k]<-mean(thetai.post[,k][p:length(thetai.post[,1])])
  }

  ## intervalle de credibilité
  
  thetha.tri <- matrix(nrow=n,ncol=Niter+1)
  
  for (k in (1:n) )
  {
    thetha.tri[k,]<-sort(thetai.post[,k])
  }

  interval.cred <- matrix(nrow=n,ncol=2)
    
 for (k in (1:n) )
  {
    interval.cred[k,]<-c(thetha.tri[k,25],thetha.tri[k,975])
  }
  
  ## Comparaison résultat expérementale et théorique

  erreur<-var(theta.exper-moyen.theta.post[1,])
  moyen.Sort.Expr<-sort(theta.exper)
  moyen.Sort.Post<-sort(moyen.theta.post[1,])
  plot(moyen.Sort.Expr,moyen.Sort.Post)
  
  