#################################################
###           GRADIENT FUNCTIONS              ###
#################################################

sourceCpp("cpp_gradients.cpp")

grad_ab <- function(Phi=Phi, X, Y, Alpha, Gamma, a, b, sigma){
  
  logPi = (Phi(Alpha, Gamma, X, Y)-t(matrix(a,n,n,byrow=T))-matrix(b,n,n,byrow=T))/sigma 
  Pi = exp(logPi)
  Pitilde = Pi
  Pitilde[1,] <- rep(0,n)
  
  baseline1 = Pi*X[,1]%*% t(Y)
  baseline2 = Pi*X[,2]%*% t(Y)
  
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

grad_w <- function(Phi=Phi, X, Y, Alpha, Gamma, a, b, sigma){
  
  gradab = cpp_gradient_ab(X, Y, sigma, Alpha, Gamma, a, b)
  grada = gradab$Da
  gradb = gradab$Db
  
  argmin = OLS_wage(Phi, X, Y, Alpha, Gamma, a, b, sigma)
  
  grad_wAlpha1 = (argmin[[2]]*(grada[,1]-X[,1]*t(Y))-argmin[[1]]*gradb[,1])/sigma
  grad_wAlpha2 = (argmin[[2]]*(grada[,2]-X[,2]*t(Y))-argmin[[1]]*gradb[,2])/sigma
  
  grad_wGamma1 = (argmin[[1]]*(-gradb[,1]+X[,1]*t(Y))+argmin[[2]]*grada[,1])/sigma
  grad_wGamma2 = (argmin[[1]]*(-gradb[,2]+X[,2]*t(Y))+argmin[[2]]*grada[,2])/sigma
  
  output <- list(grad_wAlpha1, grad_wAlpha2, grad_wGamma1, grad_wGamma2)
  return(output)
}

grad_L2 <- function(Phi=Phi, X, Y, Alpha, Gamma, a, b, sigma){
  
  gradw = grad_w(Phi, X, Y, Alpha, Gamma, a, b, sigma)
  
  argmin = OLS_wage(Phi, X, Y, Alpha, Gamma, a, b, sigma)
  
  res = argmin[[5]]/argmin[[4]]
  
  L2Alpha1 = res%*%t(gradw[[1]])
  L2Alpha2 = res%*%t(gradw[[2]])
  L2Gamma1 = res%*%t(gradw[[3]])
  L2Gamma2 = res%*%t(gradw[[4]])
  
  output <- list(L2Alpha1, L2Alpha2, L2Gamma1, L2Gamma2)
  return(output)
  
}

grad_L1 <- function(Phi=Phi, X, Y, Alpha, Gamma, a, b, sigma){
  
  gradab = cpp_gradient_ab(X, Y, sigma, Alpha, Gamma, a, b)
  grada = gradab$Da
  gradb = gradab$Db
  
  L1Alpha1 = sum(X[,1]*Y - grada[,1] - gradb[,1])/sigma
  L1Alpha2 = sum(X[,2]*Y - grada[,2] - gradb[,2])/sigma
  
  output <- list(L1Alpha1, L1Alpha2, L1Alpha1, L1Alpha2)
  return(output)
}


# grad_L1 <- function(Phi=Phi, X, Y, Alpha, Gamma, sigma = 1){
#   
#   output_IPFP = IPFP(Phi, X, Y, Alpha, Gamma, sigma)
#   a = output_IPFP[[1]]
#   b = output_IPFP[[2]]
#   
#   logPi = (Phi(Alpha, Gamma, X, Y)-matrix(a,n,n,byrow=T)-t(matrix(b,n,n,byrow=T)))/sigma 
#   Pi = exp(logPi)
#   
#   L1Alpha1 = (sum(diag(X[,1]%*%t(Y))) - sum(Pi*X[,1]%*%t(Y)))/sigma
#   L1Alpha2 = (sum(diag(X[,2]%*%t(Y))) - sum(Pi*X[,2]%*%t(Y)))/sigma
#   
#   output <- list(L1Alpha1, L1Alpha2, L1Alpha1, L1Alpha2)
#   return(output)
# }

grad_allin1 <- function(x){
  
  Alpha <- x[1:2]
  Gamma <- x[3:4]
  
  outputIPFP = cpp_IPFP_lse(Xvals, Yvals, sigma, Alpha, Gamma, 1E4, 1E-8)
  
  a = outputIPFP$a
  b = outputIPFP$b
  
  output = unlist(grad_L1(Phi, Xvals, Yvals, Alpha, Gamma, a, b, sigma)) + unlist(grad_L2(Phi, Xvals, Yvals, Alpha, Gamma, a, b, sigma))
  return(output)
}
