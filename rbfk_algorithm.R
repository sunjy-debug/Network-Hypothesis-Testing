.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")
library(noisySBM)

rbfk_algorithm = function(file_directory, output_directory, data_generation_method, random_seed, tau = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.25)){
  data = readRDS(file_directory)
  X = data$X
  A = data$A
  Z = data$Z
  
  set.seed(random_seed)
  
  start_rbfk = Sys.time()
  cat("Rebafka's algorithm starts.\n")
  fit = noisySBM::fitNSBM(X, model = "Gauss01", initParam = list(maxNbOfPasses = 1))
  infer_rbfk = lapply(tau, function(x) noisySBM::graphInference(X, nodeClustering = fit_rbfk[[2]]$clustering, theta = fit_rbfk[[2]]$theta, alpha = x, modelFamily = "Gauss")$A)
  TDR_rbfk =  sapply(infer_rbfk, function(x) sum(A * x, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE))
  FDR_rbfk = sapply(infer_rbfk, function(x) sum((1 - A) * x, na.rm = TRUE) / pmax(sum(x, na.rm = TRUE), 1))
  rindex_rbfk = adj.rand.index(fit_rbfk[[2]]$clustering, Z)
  end_rbfk = Sys.time()
  elapse_rbfk = end_rbfk - start_rbfk
  cat("Rebafka's algorithm ends: \n")
  print(elapse_rbfk)
  
  return(list(TDR_rbfk = TDR_rbfk, FDR_rbfk = FDR_rbfk, TIME_rbfk = elapse_rbfk))
}

rbfk_main = function(){
  library(rslurm)
  indir <- Sys.getenv("INDIR", unset = file.path(getwd(), "Data"))
  outdir = Sys.getenv("OUTDIR", unset = file.path(getwd(), "Results/RBFK"))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  options(bitmapType = "cairo")
  
  niteration = 500
  cat("Submitting Slurm job...\n")
  sjob = slurm_apply(
    rbfk_algorithm,
    params = data.frame(file_directory = file.path(indir, paste0("NSBM_", 1:niteration, ".rds")), output_directory = outdir, data_generation_method = "NSBM", random_seed = 1:niteration),
    jobname = "rbfk",
    nodes = niteration,
    cpus_per_node = 1,
    pkgs = c(
      "noisySBM"
    )
  )
  cat("Slurm job submitted!\n")
  results = get_slurm_out(sjob, outtype = "raw")
}