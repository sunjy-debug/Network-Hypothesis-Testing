.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")
library(igraph)

## data generation
# NSBM setting
NSBM_generation = function(){
  n = 20
  k = 2
  Z = sample.int(k, size = n, replace = TRUE, prob = rep(1 / 2, k))
  Q = matrix(0.1, k, k)
  diag(Q) = 0.8
  mu = matrix(3, k, k)
  diag(mu) = 1
  sigma = matrix(1, k, k)
  A = matrix(0, n, n)
  for (i in 1:n){
    r = Z[i]
    for (j in i:n){
      s = Z[j]
      if(i == j) next
      A[i,j] = rbinom(1, 1, prob = Q[r,s])
      A[j,i] = A[i,j]
    }
  }
  X = matrix(0, n, n)
  for (i in 1:n){
    r = Z[i]
    for (j in i:n){
      s = Z[j]
      if(A[i,j] == 1){
        X[i,j] = rnorm(1, mu[r,s], sigma[r,s])
      } else {
        X[i,j] = rnorm(1, 0, 1)
      }
      X[j,i] = X[i,j]
    }
  }
  return(list(A = A, X = X, Z = Z))
}

# star A setting
stargraph = function(){
  n = 200
  A = matrix(0, n, n)
  A[1, -1] = 1
  A[-1, 1] = 1
  X = matrix(0, n, n)
  for (i in 1:n){
    for (j in i:n){
      if(A[i,j] == 1){
        X[i,j] = rnorm(1, 2, 1)
      } else {
        X[i,j] = rnorm(1, 0, 1)
      }
      X[j,i] = X[i,j]
    }
  }
  return(list(A = A, X = X))
}

# random bipartite A setting
randombipartitegraph = function(){
  n = 200
  group1 = sample.int(n, size = floor(n / 2), replace = FALSE)
  group2 = setdiff(1:n, group1)
  A = matrix(0, n, n)
  for(i in group1){
    for(j in group2){
      A[i, j] = rbinom(1, 1, .5)
      A[j, i] = A[i, j]
    }
  }
  X = matrix(0, n, n)
  for (i in 1:n){
    for (j in i:n){
      if(A[i,j] == 1){
        X[i,j] = rnorm(1, 2, 1)
      } else {
        X[i,j] = rnorm(1, 0, 1)
      }
      X[j,i] = X[i,j]
    }
  }
  return(list(A = A, X = X))
}

# preferential attachment A setting
preferattachgraph = function(){
  n = 200
  root_graph = erdos.renyi.game(n = 8, p = 0.5, type = "gnp") # root graph (Erdos-Renyi graph)
  pa_graph = sample_pa(n = n, power = 1, m = 6, out.pref = TRUE, start.graph = root_graph, directed = FALSE)
  A = as_adjacency_matrix(pa_graph, type = "both", sparse = FALSE)
  
  X = matrix(0, n, n)
  for (i in 1:n){
    for (j in i:n){
      if(A[i,j] == 1){
        X[i,j] = rnorm(1, 2, 1)
      } else {
        X[i,j] = rnorm(1, 0, 1)
      }
      X[j,i] = X[i,j]
    }
  }
  return(list(A = A, X = X))
}

data_generation = function(output_directory, data_generation_method, random_seed){
  set.seed(random_seed)
  
  start_data = Sys.time()
  cat("Data generation starts.\n")
  if(data_generation_method == "NSBM"){
    data_generation = NSBM_generation()
  }else if (data_generation_method == "star"){
    data_generation = stargraph()
  }else if (data_generation_method == "randombi"){
    data_generation = randombipartitegraph()
  }else if (data_generation_method == "preferattach"){
    data_generation = preferattachgraph()
  }
  X = data_generation$X
  A = data_generation$A
  if(data_generation_method == "NSBM"){
    Z = data_generation$Z
  } else {
    Z = NULL
  }
  end_data = Sys.time()
  elapse_data = end_data - start_data
  cat("Data generation ends: \n")
  print(elapse_data)
  
  if (!dir.exists(output_directory)) dir.create(output_directory, recursive = TRUE, showWarnings = FALSE)
  data = list(X = X, A = A, Z = Z)
  path = file.path(output_directory, paste0(paste(c(data_generation_method, random_seed), collapse = "_"), ".rds"))
  saveRDS(data, path)
  cat("Data is saved. \n")
  return(invisible(NULL))
}

data_generation_main = function(){
  library(rslurm)
  outdir = Sys.getenv("OUTDIR", unset = file.path(getwd(), "Data"))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  options(bitmapType = "cairo")
  
  niteration = 500
  cat("Submitting Slurm job...\n")
  sjob = slurm_apply(
    data_generation,
    params = data.frame(output_directory = outdir, data_generation_method = "NSBM", random_seed = 1:niteration),
    jobname = "data_generation",
    nodes = niteration,
    cpus_per_node = 1,
    global_objects = c(
      "NSBM_generation", "stargraph", "randombipartitegraph", "preferattachgraph"
    ),
    pkgs = c("igraph")
  )
  cat("Slurm job submitted!\n")
  results = get_slurm_out(sjob, outtype = "raw")
}