.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")
library(noisysbmGGM)
env_noisysbmGGM = asNamespace("noisysbmGGM")
env_klln = new.env(parent = env_noisysbmGGM)
sys.source("noisysbmGGM_decouple.R", envir = env_klln)
main_noisySBM_fit = get("main_noisySBM_fit",   envir = env_klln)
main_noisySBM_infer = get("main_noisySBM_infer", envir = env_klln)
environment(main_noisySBM_fit) = env_klln
environment(main_noisySBM_infer) = env_klln

klln_algorithm = function(file_directory, output_directory, data_generation_method, random_seed, tau = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.25)){
  data = readRDS(file_directory)
  X = data$X
  A = data$A
  Z = data$Z
  
  set.seed(random_seed)
  
  start_klln = Sys.time()
  cat("Kilian's algorithm starts.\n")
  fit_klln = main_noisySBM_fit(X, NIG = TRUE)
  infer_klln = lapply(tau, function(x) main_noisySBM_infer(fit_klln$Rho, fit_klln$theta, fit_klln$Z, X, alpha = x)$A)
  TDR_klln =  sapply(infer_klln, function(x) sum(A * x, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE))
  FDR_klln = sapply(infer_klln, function(x) sum((1 - A) * x, na.rm = TRUE) / pmax(sum(x, na.rm = TRUE), 1))
  rindex_klln = adj.rand.index(fit_klln$Z, Z)
  end_klln = Sys.time()
  elapse_klln = end_klln - start_klln
  cat("Killian's algorithm ends: \n")
  print(elapse_klln)
  
  return(list(TDR_klln = TDR_klln, FDR_klln = FDR_klln, TIME_klln = elapse_klln))
}

klln_main = function(){
  library(rslurm)
  indir <- Sys.getenv("INDIR", unset = file.path(getwd(), "Data"))
  outdir = Sys.getenv("OUTDIR", unset = file.path(getwd(), "Results/KLLN"))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  options(bitmapType = "cairo")
  
  niteration = 500
  cat("Submitting Slurm job...\n")
  sjob = slurm_apply(
    klln_algorithm,
    params = data.frame(file_directory = file.path(indir, paste0("NSBM_", 1:niteration, ".rds")), output_directory = outdir, data_generation_method = "NSBM", random_seed = 1:niteration),
    jobname = "klln",
    nodes = niteration,
    cpus_per_node = 1,
    global_objects = c(
      "main_noisySBM_fit", "main_noisySBM_infer"
    ),
    pkgs = c(
      "noisysbmGGM"
    )
  )
  cat("Slurm job submitted!\n")
  results = get_slurm_out(sjob, outtype = "raw")
}