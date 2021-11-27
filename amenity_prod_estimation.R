#######################################################
############ Pauline Corblet - Sciences Po ############
#######################################################
#
#######################################################
### Job Amenities and Labor Productivity Estimation ###
############## From Dupuy & Galichon 2017 #############
#######################################################

#######################################################
##################### DATA PREP #######################
#######################################################

library(MASS)
#library(Brobdingnag)
#library(gradDescent)
library(Matrix)
library(gurobi)
library(neldermead)

seed = 777
set.seed(seed)

n = 5
mydata = mvrnorm(n=n, mu=c(0,0,0), Sigma = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1), nrow=3), empirical = TRUE)
#wage = runif(n=n, min = 500, max = 2000)
match = data.frame(mydata)
names(match) <-  c("employee1", "employee2","job1")

Xvals = as.matrix(data.frame(match$employee1, match$employee2))
Yvals = as.matrix(data.frame(match$job1))

###### data preprocessing
meanX = apply(Xvals,2,mean)
meanY = apply(Yvals,2,mean)
sdX = apply(Xvals,2,sd)
sdY = apply(Yvals,2,sd)

Xvals = t( t(Xvals) - meanX )
Yvals = t( t(Yvals) - meanY )
Xvals = t( t(Xvals) / sdX)
Yvals = t( t(Yvals) / sdY)

dX=dim(Xvals)[2]
dY=dim(Yvals)[2]

coeffs_Alpha = matrix(rep(1,dX*dY),nrow=dX)
coeffs_Gamma = matrix(rep(1,dX*dY),nrow=dX)
coeffs_Phi =  coeffs_Alpha+coeffs_Gamma

Phi <- function(Alpha, Gamma, X, Y) {
  dim = dim(X)[2]
  Phi = Alpha + Gamma
  return(X %*% matrix(Phi,nrow=dim) %*% t(Y))
}

##############################################################
################# OPTIMAL ASSIGNEMENT #########################
##############################################################

Phi_vec = Phi(c(1,1), c(1,1), X=Xvals, Y=Yvals)
p = rep(1,n)
q = rep(1,n)

N = dim(Phi_vec)[1]
M = dim(Phi_vec)[2]

c = c(Phi_vec)

A1 = kronecker(matrix(1,1,M),sparseMatrix(1:N,1:N))
A2 = kronecker(sparseMatrix(1:M,1:M),matrix(1,1,N))
A = rbind2(A1,A2)

d = c(p,q) 

result = gurobi(list(A=A,obj=c,modelsense="max",rhs=d,sense="="), params=list(OutputFlag=0)) 
if (result$status=="OPTIMAL") {
  pi = matrix(result$x,nrow=N)
  u = result$pi[1:N]
  v = result$pi[(N+1):(N+M)]
  val = result$objval
} else {stop("optimization problem with Gurobi.") }


# Set up otimal matches
match = which(pi==1, arr.ind=TRUE)
match = match[order(match[,1]),]

Xvals_match = Xvals
Yvals_match = as.matrix(Yvals[match[,2]])

epsilon = rnorm(n)
wage = rowSums(pi*Phi_vec)*0.5 + epsilon
#wage = rowSums(pi*Phi_vec)*0.5

##############################################################
########################## IPFP ##############################
##############################################################

IPFP <- function(Phi=Phi, Alpha, Gamma, sigma, tolIpfp=1E-12, maxiterIpfp=1E6, B = rep(1,n)){

  f = rep(1/n,n)
  g = rep(1/n,n)

  #B=rep(1,n)
  C=exp(Phi(Alpha, Gamma, Xvals_match, Yvals_match)/sigma)

  contIpfp = TRUE
  iterIpfp = 0

  while(contIpfp)
  {
    iterIpfp = iterIpfp+1

    A = f / (C %*% B)
    CA = c(t(A) %*% C)
    B = matrix(g / CA, nrow =n)

    if (((max(abs(CA*B - g))<tolIpfp ) & (max(abs(A*(c(C%*%B)) - f))<tolIpfp )) | (iterIpfp >= maxiterIpfp)) {contIpfp=FALSE}

  }

  norm = A[1]
  A = A/norm
  B = B*norm

  u = -sigma*log(A)
  v = -sigma*log(B)

  output <- list(u,v,iterIpfp)
  return(output)
}


IPFP <- function(Phi=Phi, Alpha, Gamma, sigma, tolIpfp=1E-17, maxiterIpfp=1E6){

  f = rep(1/n,n)
  g = rep(1/n,n)

  v = rep(0,n)
  mu = - sigma * log(f)
  nu = - sigma * log(g)
  uprec = -Inf

  surplus = Phi(Alpha, Gamma, Xvals_match, Yvals_match)

  contIpfp = TRUE
  iterIpfp = 0

  while(contIpfp)
  {
    iterIpfp = iterIpfp+1

    vstar = apply(surplus - matrix(v,n,n,byrow=T), 1, max)

    u=mu + vstar + sigma * log( apply(exp( (surplus - t(matrix(v,n,n,byrow=T)) - matrix(vstar,n,n,byrow=T))/sigma) ,1, sum) )
    error = max(abs(u-uprec))
    uprec = u

    ustar = apply(surplus - u, 2, max)
    v = nu + ustar + sigma * log(apply(exp( (surplus - matrix(u,n,n,byrow=T) - t(matrix(ustar,n,n,byrow=T)) )/sigma) ,2, sum))

    #if( (error<tolIpfp) | (iterIpfp >= maxiterIpfp)) {contIpfp=FALSE}
    test = surplus-matrix(u,n,n,byrow=T)-t(matrix(v,n,n,byrow=T))
    if( ( max(abs(rowSums(test)-f)) < tolIpfp & max(abs(colSums(test)-g)) < tolIpfp)  | (iterIpfp >= maxiterIpfp)) {contIpfp=FALSE}
  }

  A = exp(-u)
  B = exp(-v)

  norm = A[1]
  A = A/norm
  B = B*norm

  u = -log(A)
  v = -log(B)

  output <- list(u,v,iterIpfp)
  return(output)
}

a = IPFP(Phi, coeffs_Alpha, coeffs_Gamma, sigma = 1)[[1]]
b = IPFP(Phi, coeffs_Alpha, coeffs_Gamma, sigma = 1)[[2]]

rowSums(exp((Phi(coeffs_Alpha, coeffs_Gamma, Xvals_match, Yvals_match)-t(matrix(a,n,n,byrow=T))-matrix(b,n,n,byrow=T)))) 
colSums(exp((Phi(coeffs_Alpha, coeffs_Gamma, Xvals_match, Yvals_match)-t(matrix(a,n,n,byrow=T))-matrix(b,n,n,byrow=T))))  


##############################################################
######################## OLS WAGE ############################
##############################################################

OLS_wage <- function(Phi=Phi, Alpha, Gamma, a, b, sigma){
  
  gamma_diag = (diag(Phi(matrix(rep(0,dX*dY),nrow=dX), Gamma, Xvals_match, Yvals_match)) - b)/sigma
  alpha_diag = (a - diag(Phi(Alpha, matrix(rep(0,dX*dY),nrow=dX), Xvals_match, Yvals_match)))/sigma
  
  OLS_wages = lm(wage ~ gamma_diag + alpha_diag)
  summary(OLS_wages)
  coeff = OLS_wages$coefficients
  
  sigma1=coeff[2]
  sigma2=coeff[3]
  t=coeff[1]
  
  ssquare = sum(OLS_wages$residuals^2)/n 

  output <- list(sigma1, sigma2, t, ssquare, OLS_wages$residuals)
  
  return(output)
  
}

##############################################################
###################### GRADIENTS SET UP ######################
##############################################################

grad_ab <- function(Phi=Phi, Alpha, Gamma, a, b, sigma){
  
  logPi = (Phi(Alpha, Gamma, Xvals_match, Yvals_match)-t(matrix(a,n,n,byrow=T))-matrix(b,n,n,byrow=T))/sigma 
  Pi = exp(logPi)
  Pitilde = Pi
  Pitilde[1,] <- rep(0,n)
  
  baseline1 = Pi*Xvals_match[,1]%*% t(Yvals_match)
  baseline2 = Pi*Xvals_match[,2]%*% t(Yvals_match)
  
  E = cbind(matrix(rowSums(baseline1), nrow = n),matrix(rowSums(baseline2), nrow = n))
  E[1,] <- 0
  F = cbind(matrix(colSums(baseline1), nrow= n),matrix(colSums(baseline2), nrow= n))
  
  block = rbind(cbind(diag(1/n,n),Pitilde),cbind(t(Pi),diag(1/n,n)))
  
  DaDb = solve(block)%*%rbind(E,F)

  m=n+1
  l=2*n
  output <- list(DaDb[1:n,], DaDb[m:l,])
  
  return(output)
}

grad_w <- function(Phi=Phi, Alpha, Gamma, a, b, sigma){
  
  gradab = grad_ab(Phi, Alpha, Gamma, a, b, sigma)
  grada = gradab[[1]]
  gradb = gradab[[2]]
  
  argmin = OLS_wage(Phi, Alpha, Gamma, a, b, sigma)
  
  grad_wAlpha1 = (argmin[[2]]*(grada[,1]-Xvals_match[,1]*t(Yvals_match))-argmin[[1]]*gradb[,1])/sigma
  grad_wAlpha2 = (argmin[[2]]*(grada[,2]-Xvals_match[,2]*t(Yvals_match))-argmin[[1]]*gradb[,2])/sigma
  
  grad_wGamma1 = (argmin[[1]]*(-gradb[,1]+Xvals_match[,1]*t(Yvals_match))+argmin[[2]]*grada[,1])/sigma
  grad_wGamma2 = (argmin[[1]]*(-gradb[,2]+Xvals_match[,2]*t(Yvals_match))+argmin[[2]]*grada[,2])/sigma
  
  output <- list(grad_wAlpha1, grad_wAlpha2, grad_wGamma1, grad_wGamma2)
  return(output)
}

grad_L2 <- function(Phi=Phi, Alpha, Gamma, a, b, sigma){
  
  gradw = grad_w(Phi, Alpha, Gamma, a, b, sigma)

  argmin = OLS_wage(Phi, Alpha, Gamma, a, b, sigma)
  
  res = argmin[[5]]/argmin[[4]]
  
  L2Alpha1 = res%*%t(gradw[[1]])
  L2Alpha2 = res%*%t(gradw[[2]])
  L2Gamma1 = res%*%t(gradw[[3]])
  L2Gamma2 = res%*%t(gradw[[4]])
  
  output <- list(L2Alpha1, L2Alpha2, L2Gamma1, L2Gamma2)
  return(output)
  
}

grad_L1 <- function(Phi=Phi, Alpha, Gamma, a, b, sigma){
  
  gradab = grad_ab(Phi, Alpha, Gamma, a, b, sigma)
  grada = gradab[[1]]
  gradb = gradab[[2]]
  
  L1Alpha1 = sum(Xvals_match[,1]*Yvals_match - grada[,1] - gradb[,1])/sigma
  L1Alpha2 = sum(Xvals_match[,2]*Yvals_match - grada[,2] - gradb[,2])/sigma
  
  output <- list(L1Alpha1, L1Alpha2, L1Alpha1, L1Alpha2)
  return(output)
}


# grad_L1 <- function(Phi=Phi, Alpha, Gamma, sigma = 1){
#   
#   output_IPFP = IPFP(Phi, Alpha, Gamma, sigma)
#   a = output_IPFP[[1]]
#   b = output_IPFP[[2]]
#   
#   logPi = (Phi(Alpha, Gamma, Xvals_match, Yvals_match)-matrix(a,n,n,byrow=T)-t(matrix(b,n,n,byrow=T)))/sigma 
#   Pi = exp(logPi)
#   
#   L1Alpha1 = (sum(diag(Xvals_match[,1]%*%t(Yvals_match))) - sum(Pi*Xvals_match[,1]%*%t(Yvals_match)))/sigma
#   L1Alpha2 = (sum(diag(Xvals_match[,2]%*%t(Yvals_match))) - sum(Pi*Xvals_match[,2]%*%t(Yvals_match)))/sigma
#   
#   output <- list(L1Alpha1, L1Alpha2, L1Alpha1, L1Alpha2)
#   return(output)
# }



##############################################################
################# LOG LIKELIHOOD MAXIMIZATION ################
##############################################################

Log_likelihood <- function(Alpha, Gamma, a, b, sigma){
  
  argmin = OLS_wage(Phi, Alpha, Gamma, a, b, sigma)

  L1 = sum(diag(Phi(Alpha,Gamma, Xvals_match, Yvals_match))-a-b)/sigma
  L2 = -sum(argmin[[5]]^2)/(2*argmin[[4]])-n*log(argmin[[4]])/2
  
  output = c(L1+L2, L1, L2)
  return(output)
}


Log_likelihood_allinc2 <- function(x){
  
  output_IPFP = IPFP(Phi, c(x,y), Gamma, sigma)
  a = output_IPFP[[1]]
  b = output_IPFP[[2]]
  
  argmin = OLS_wage(Phi, c(x,y), Gamma, a, b, sigma)
  
  L1 = sum(diag(Phi(c(x,y),Gamma, Xvals_match, Yvals_match))-a-b)/sigma
  L2 = -sum(argmin[[5]]^2)/(2*argmin[[4]])-n*log(argmin[[4]])/2
  
  return(L2)
}

tolGD = 1E-14
maxiterGD = 1E6
stepGD = 1E-2

contGD = TRUE
iterGD = 0

Alpha = c(13.3,1)
Gamma = c(1,1)

scale = 1

B_prev = rep(1,n)

while(contGD){
  
  #print(Sys.time())
  iterGD = iterGD+1
  
  output_IPFP = IPFP(Phi, Alpha, Gamma, scale, maxiterIpfp = 1E7, B = B_prev)
  if (output_IPFP[[3]] == 1E7 ) {stop('maximum number of IPFP iterations reached')}
  a = output_IPFP[[1]]
  b = output_IPFP[[2]]
  
  B_prev = exp(-b/sigma)
  
  #print(paste0("Current L1: ",Log_likelihood(Alpha, Gamma, a, b, scale)[2]))
  print(paste0("Current L2: ",Log_likelihood(Alpha, Gamma, a, b, sigma = scale)[3]))
  #print(paste0("Current L1+L2: ",Log_likelihood(Alpha, Gamma, a, b, sigma = scale)[1]))
  
  thegrad1 = grad_L1(Phi,Alpha,Gamma, a, b, scale) 
  #thegrad2 = grad_L2(Phi,Alpha,Gamma, a, b, scale)
  
  epsilon = 1E-3
  target = Alpha[1]
  thegrad2 = list(Log_likelihood_allinc2((target+epsilon)-Log_likelihood_allinc2(target-epsilon))/2*epsilon,0,0,0)
  
  grad_diff <- vector("list",4)
  for (i in 1:4){
    grad_diff[[i]] <- -(thegrad2[[i]])
    #grad_diff[[i]] <- -(thegrad1[[i]]+thegrad2[[i]])
  }
  
  spread = max(abs(unlist(grad_diff)))
  #spread = abs(unlist(grad_diff)[1])
  
  print(paste0("Current number of iterations: ", iterGD))
  print(paste0(c("Current Alpha: ",Alpha), collapse= " "))
  print(paste0(c("Current Gamma: ",Gamma), collapse= " "))
  print(paste0(c("Current gradient: ",grad_diff),collapse= " "))
  print(" ")
  
  if (iterGD >= maxiterGD ) {stop('maximum number of iterations reached')}
  if (spread < tolGD | iterGD >= maxiterGD ) {contGD=FALSE}
  
  #Alpha = Alpha - stepGD*unlist(grad_diff)[1:2]
  #Alpha[2] = Alpha[2]- stepGD*unlist(grad_diff)[2]
  #Gamma = Gamma - stepGD*unlist(grad_diff)[3:4]
  Alpha[1] = Alpha[1]- stepGD*unlist(grad_diff)[1]
}


##############################################################
############################ PLOTS ###########################
##############################################################

start = c(10,1)
Gamma = c(1,1)
scale = 10

output_IPFP = IPFP(Phi, start, Gamma, scale)
a = output_IPFP[[1]]
b = output_IPFP[[2]]

values = as.array(Log_likelihood(start, Gamma, a,b,sigma = scale)[[1]])

for(i in seq(from=10.01, to=20, by=0.01)){
  Alpha = matrix(c(i,1), nrow=dX)
  
  output_IPFP = IPFP(Phi, Alpha, Gamma, scale)
  a = output_IPFP[[1]]
  b = output_IPFP[[2]]
  
  values = cbind(values,Log_likelihood(Alpha, Gamma, a, b, sigma =scale)[[1]])
}

plot(seq(from=10, to=20, by=0.01),values, xlab="Parameter value", ylab="Log likelihood", main="",pch=20)



y = 1
Gamma = c(1,1)
sigma = 1

Log_likelihood_allinc <- function(x){
  
  output_IPFP = IPFP(Phi, c(x,y), Gamma, sigma)
  a = output_IPFP[[1]]
  b = output_IPFP[[2]]
  
  argmin = OLS_wage(Phi, c(x,y), Gamma, a, b, sigma)
  
  L1 = sum(diag(Phi(c(x,y),Gamma, Xvals_match, Yvals_match))-a-b)/sigma
  L2 = -sum(argmin[[5]]^2)/(2*argmin[[4]])-n*log(argmin[[4]])/2
  
  return(L2)
}



res=optimize(Log_likelihood_allinc, lower = 0, upper =100, maximum = TRUE)

output_IPFP = IPFP(Phi, c(res$maximum,y), Gamma, sigma)
a = output_IPFP[[1]]
b = output_IPFP[[2]]

unlist(grad_L1(Phi,c(res$maximum,y),Gamma, a, b, sigma))+unlist(grad_L2(Phi,c(res$maximum,y),Gamma, a, b, sigma)) 

values <- array()
for(i in seq(from=1, to=2, by=0.001)){
  Alpha = c(i,1)
  values = cbind(values,grad_L1(Phi,Alpha, Gamma, sigma)[[1]])
}
plot(seq(from=1, to=2, by=0.001),values[-1], xlab="Parameter value", ylab="LogL1 gradient", main="",pch=20)

values <- array()
for(i in seq(from=0, to=100, by=0.1)){
  values = cbind(values,Log_likelihood_allinc(i))
}
plot(seq(from=0, to=100, by=0.1),values[-1], xlab="Parameter value", ylab="Log likelihood", main="",pch=20)



sigma =1
Log_likelihood_allin1 <- function(x){
  
  Alpha <- x[1:2]
  Gamma <- x[3:4]
  
  output_IPFP = IPFP(Phi, Alpha, Gamma, sigma)
  a = output_IPFP[[1]]
  b = output_IPFP[[2]]
  
  argmin = OLS_wage(Phi, Alpha, Gamma, a, b, sigma)
  
  L1 = sum(diag(Phi(Alpha,Gamma, Xvals_match, Yvals_match))-a-b)/sigma
  L2 = -sum(argmin[[5]]^2)/(2*argmin[[4]])-n*log(argmin[[4]])/2
  
  return(L1+L2)
}


res = optim(par=c(1,1,1,1),fn=Log_likelihood_allin1, control = list(trace = 1, fnscale = -1, maxit = 1E3))
Alpha = res$par[1:2]
Gamma = res$par[3:4]