Niter=5
n=12
thetai.post <- matrix(nrow=Niter,ncol=n)
for(j in 1:Niter )
{  
for(k in 1:n )
{
  # beta(yi+a,ni-yi+b)
  thetai.post[j,k]=rbeta(1,shape1=1,shape2=2)
}
}