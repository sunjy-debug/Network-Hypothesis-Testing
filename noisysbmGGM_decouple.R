rm(list=ls())

.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")
library(noisysbmGGM)

main_noisySBM_fit = function (X, NIG = FALSE, threshold = 0.5, Nbrepet = 2, rho = NULL, 
                              tau = NULL, a = NULL, b = NULL, c = NULL, d = NULL, n0 = 1, 
                              eta0 = 1, zeta0 = 1, Qup = NULL, nbCores = parallel::detectCores(), 
                              nbOfZ = 12, sigma0 = 1, sigma1 = 1, percentageOfPerturbation = 0.3, 
                              verbatim = TRUE) 
{
  p = length(X[1, ])
  if (is.null(Qup)) 
    Qup = 10
  if (length(X[, 1]) != p) {
    stop("X must be a square matrix")
  }
  if (NIG) {
    if (is.null(a)) 
      a = 0
    if (is.null(b)) 
      b = 1
    if (is.null(c)) 
      c = 1
    if (is.null(d)) 
      d = 1
    dataVec = X[lower.tri(X)]
    result = mainSearchOpti_NIG(X, Qup = Qup, threshold = threshold, 
                                Nbrepet = Nbrepet, nbCores = nbCores, nbOfZ = nbOfZ, 
                                percentageOfPerturbation = percentageOfPerturbation, 
                                a = a, b = b, c = c, d = d, n0 = n0, eta0 = eta0, 
                                zeta0 = zeta0, sigma0 = sigma0, fast = TRUE, verbatim = verbatim)
    if (result$Q > 1) {
      test = Merge_Nondirige_NIG(dataVec, result$Z, result$Zmatrix, 
                                 result$Q, result$Rho, result$A, result$theta, 
                                 threshold = threshold, a = a, b = b, c = c, 
                                 d = d, n0 = n0, eta0 = eta0, zeta0 = zeta0)
      return(list(Rho = test$Rhomerge, theta = test$thetamerge, Z = test$Zmerge))
      
    }
    else {
      return(list(Rho = result$Rho, theta = result$theta, Z = result$Z))
    }
  }
  else {
    if (is.null(rho)) 
      rho = 1
    if (is.null(tau)) 
      tau = 1
    dataVec = X[lower.tri(X)]
    result = mainSearchOpti(X, Qup = Qup, threshold = threshold, 
                            Nbrepet = Nbrepet, nbCores = nbCores, nbOfZ = nbOfZ, 
                            percentageOfPerturbation = percentageOfPerturbation, 
                            rho = rho, tau = tau, n0 = n0, eta0 = eta0, zeta0 = zeta0, 
                            sigma0 = sigma0, sigma1 = sigma1, fast = TRUE, verbatim = verbatim)
    if (result$Q > 1) {
      test = Merge_Nondirige(dataVec, result$Z, result$Zmatrix, 
                             result$Q, result$Rho, result$A, result$theta, 
                             threshold = threshold, rho = rho, tau = tau, 
                             n0 = n0, eta0 = eta0, zeta0 = zeta0)
      return(list(Rho = test$Rhomerge, theta = test$thetamerge, Z = test$Zmerge))
    }
    else {
      return(list(Rho = result$Rho, theta = result$theta, Z = result$Z))
    }
  }
}

main_noisySBM_infer = function (Rho, theta, Z, X, alpha) 
{
  dataVec = X[lower.tri(X)]
  p = length(Z)
  A = matrix(0, p, p)
  Avec = get_Aall_qval(Rho, theta, Z, dataVec, alpha)
  A[lower.tri(A)] = Avec
  A <- A + t(A)
  return(A = A)
}