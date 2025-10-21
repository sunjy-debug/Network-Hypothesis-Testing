.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")
library(microbenchmark)
library(noisySBM)
library(noisysbmGGM)
source("noisysbmGGM_decouple.R")
env_mfmsbmseq = new.env(parent = emptyenv())
env_mfmsbmsim = new.env(parent = emptyenv())
env_sgdpmmsbmseq = new.env(parent = emptyenv())
env_sgdpmmsbmsim = new.env(parent = emptyenv())

source("MFMSBM_sequential.R", local = env_mfmsbmseq)
source("MFMSBM_simultaneous.R", local = env_mfmsbmsim)
source("SGDPMMSBM_sequential.R", local = env_sgdpmmsbmseq)
source("SGDPMMSBM_simultaneous.R", local = env_sgdpmmsbmsim)

set.seed(0)
data_generation = env_mfmsbmseq$NSBM_generation()
X = data_generation$X
tau = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.25)
mb = microbenchmark(
  MFMSBM_seq = {
    fit_MFM = env_mfmsbmseq$MFMSBM_estimation(X, niterations = 300, delta = 2, xi = .5, rou = 2.5, kappa = 1, alpha = 4, beta = 7, gamma = 10, lambda = 2)
    infer_MFM = lapply(tau, function(x) env_mfmsbmseq$MFMSBM_inference(fit_MFM$q, x))
    invisible(infer_MFM)
  },
  MFMSBM_sim = {
    fit_MFM = env_mfmsbmsim$MFMSBM_estimation(X, niterations = 300, delta = 2, xi = .5, rou = 2.5, kappa = 1, alpha = 4, beta = 7, gamma = 10, lambda = 2)
    infer_MFM = lapply(tau, function(x) env_mfmsbmsim$MFMSBM_inference(fit_MFM$q, x))
    invisible(infer_MFM)
  },
  SGDPMMSBM_seq = {
    fit_SGDPMM = env_sgdpmmsbmseq$SGDPMMSBM_estimation(X, niterations = 300, delta = 5, xi = 2.5, rou = 1.2, kappa = 15, alpha = 4, beta = 6, eta = 2.5, zeta = 2, m = 5)
    infer_SGDPMM = lapply(tau, function(x) env_sgdpmmsbmseq$SGDPMMSBM_inference(fit_SGDPMM$q, x))
    invisible(infer_SGDPMM)
  },
  SGDPMMSBM_sim = {
    fit_SGDPMM = env_sgdpmmsbmsim$SGDPMMSBM_estimation(X, niterations = 300, delta = 5, xi = 2.5, rou = 1.2, kappa = 15, alpha = 4, beta = 6, eta = 2.5, zeta = 2, m = 5)
    infer_SGDPMM = lapply(tau, function(x) env_sgdpmmsbmsim$SGDPMMSBM_inference(fit_SGDPMM$q, x))
    invisible(infer_SGDPMM)
  },
  RBFK = {
    fit_rbfk = noisySBM::fitNSBM(X, model = "Gauss01")
    infer_rbfk = lapply(tau, function(x) noisySBM::graphInference(X, nodeClustering = fit_rbfk[[2]]$clustering, theta = fit_rbfk[[2]]$theta, alpha = x, modelFamily = "Gauss")$A)
    invisible(infer_rbfk)
  },
  KLLN = {
    fit_klln = main_noisySBM_fit(X, NIG = TRUE)
    infer_klln = lapply(tau, function(x) main_noisySBM_infer(fit_klln$Rho, fit_klln$theta, fit_klln$Z, X, alpha = x)$A)
    invisible(infer_klln)
  },
  times = 3,
  unit  = "ms"
)
results = summary(mb)
output_dir = "Output"
write.csv(results, file.path(output_dir, "results.csv"), row.names = FALSE)