.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")
library(dplyr)
library(tidyr)
library(ggplot2)
library(fossil)
library(invgamma)
library(MASS)
library(igraph)
library(noisySBM)
library(noisysbmGGM)
source("noisysbmGGM_decouple.R")

## Function for log-likelihood related to jth observation
loglike = function(Z, Q, mu, sigma, X, j, n)
{
  ################################################################
  
  ## Input: Z = clustering configuration, a n by 1 vector ##
  ##        Q = probability matrix, a k by k matrix ##
  ##        mu = mean matrix, a k by k matrix ##
  ##        sigma = standard deviation matrix, a k by k matrix ##
  ##        X = the observation data matrix, a n by n matrix ##
  ##        j = observation index ##
  ##        n = number of observations ##
  
  ## output: log-likelihood related to Jth observation ##
  
  #################################################################
  Q = as.matrix(Q)
  
  output = 0
  if (j < n) {
    idx = (j+1):n
    output = output + sum(log((1 - Q[Z[j], Z[idx]]) * dnorm(X[j, idx], 0, 1) +
                                Q[Z[j], Z[idx]] * dnorm(X[j, idx], mu[Z[j], Z[idx]], sigma[Z[j], Z[idx]])))
  }
  if (j > 1) {
    idx = 1:(j-1)
    output = output + sum(log((1 - Q[Z[idx], Z[j]]) * dnorm(X[idx, j], 0, 1) +
                                Q[Z[idx], Z[j]] * dnorm(X[idx, j], mu[Z[idx], Z[j]], sigma[Z[idx], Z[j]])))
  }
  return(output)
}

#function for getting m(Aj)
logmargs = function(Z, X, j, alpha, beta, rou, kappa, delta, xi)
{
  ################################################################
  
  ## Input: Z = clustering configuration, a n by 1 vector ##
  ##        X = the observation data matrix, a n by n matrix ##
  ##        j = observation index ##
  ##        alpha, beta = hyperparameters for the prior on elements in Q matrix in Beta distribution ##
  ##        rou, kappa = hyperparameters for the prior on elements in mu matrix in Normal distribution ##
  ##        delta, xi = hyperparameters for the prior on elements in sigma matrix in Inverse Gamma distribution ##
  
  ## Output: m(A_j) in collapsed sampler for MFM-SBM ##
  
  #################################################################
  n = nrow(X)
  index = setdiff(1:n, j)
  x = X[cbind(pmin(index, j), pmax(index, j))]
  EQ1  = alpha / (alpha + beta) # E[Q]
  EQ0 <- beta / (alpha + beta) # E[1-Q]
  s = sqrt(xi * (kappa + 1) / (kappa * delta))
  f0 = dnorm(x, 0, 1)
  f1 = dt((x - rou)/s, df = 2 * delta) / s
  result = sum(log(EQ0 * f0 + EQ1 * f1))
  return(result)
}

## function for Collapsed sampler for MFM-SBM (main algorithm)
MFMSBM_estimation = function(X, niterations, delta, xi, rou, kappa, alpha, beta, gamma, lambda)
{
  ## Model: X_{ij}|A_{ij},Z \sim (1 - A_{ij}) N(0,1) + A_{ij} N(mu_{Z_i,Z_j}, sigma^2_{Z_i,Z_j})
  ##        A_{ij}|Z,Q \sim Ber(Q_{Z_i,Z_j}) ##
  ##        Q_{rs} \sim Beta(alpha, beta), r,s = 1,...,k ##
  ##        P(z_i = j) = \pi_j, j = 1,...,k ##
  ##        \pi \sim Dirichlet_k(gamma) ##
  ##        k-1 \sim possion(lambda) ##
  ##        sigma^2_{rs} \sim Inv-Gamma(delta, xi) ##
  ##        mu_{rs}|sigma^2_{rs} \sim N(rou, sigma^2_{rs} / kappa) ##
  
  
  
  ################################################################
  
  ## Input: X = the observation data matrix, a n by n matrix ##
  ##        niterations = the total number of iterations in MFM-SBM ##
  ##        alpha, beta = hyperparameters for the prior on elements in Q matrix in Beta distribution ##
  ##        gamma = the parameter in Dirichlet distribution that controls the relative size of clusters ##
  ##        lambda = the parameter for Poisson distribution ##
  ##        rou, tao = hyperparameters for the prior on elements in mu vector in Normal distribution ##
  ##        delta, xi = hyperparameters for the prior on elements in sigma2 matrxi in Inverse Gamma distribution ##
  
  
  ## Output: 
  ##         zout = clustering configuration, a n by 1 vector##
  ##         Qout = probability matrix, a k by k matrix ##
  
  #################################################################
  n = dim(X)[1]
  
  # initialization of clustering configuration
  # Z, k, Q, mu, sigma^2
  k = rpois(1, lambda) + 1
  Z = sample.int(k, size = n, replace = TRUE) # cluster assignment
  Z = as.integer(factor(Z))
  counts_Z = table(as.factor(Z)) # cluster sizes
  k = length(counts_Z) # number of clusters
  Q = matrix(0, k, k)
  sigma = matrix(0, k, k)
  mu = matrix(0, k, k)
  for (i in 1:k){
    for (j in i:k){
      Q[i,j] = rbeta(1, alpha, beta)
      Q[j,i] = Q[i,j]
      sigma[i,j] = sqrt(rinvgamma(1, delta, xi))
      sigma[j,i] = sigma[i,j]
      mu[i,j] = rnorm(1, rou, sigma[i,j] / sqrt(kappa))
      mu[j,i] = mu[i,j]
    }
  }
  
  # History = vector("list", niterations)
  
  # Vn: a set of pre-calculated logarithmic normalization constants related to the number of clusters
  # adjusting the probability weight of the number of clusters by encoding the prior probability structure of the number of clusters
  Vn = numeric(n + 1)
  max_w = n + 1
  for (t in 1:(n + 1)) {
    r = -Inf
    for (w in t:max_w) {
      sum_log_numerator = lgamma(w + 1) - lgamma(w - t + 1)
      sum_log_denominator = lgamma(w * gamma + n) - lgamma(w * gamma)
      b = sum_log_numerator - sum_log_denominator + dpois(w - 1, lambda, log = TRUE)
      m = max(r, b)
      r = log(exp(r - m) + exp(b - m)) + m
    }
    Vn[t] = r
  }
  
  ##start Gibb's sampling
  for (niter in 1:niterations)
  {
    Z_old = Z
    counts_Z_old = table(as.factor(Z_old))
    k_old = length(counts_Z_old)
    Q_old = Q
    mu_old = mu
    sigma_old = sigma
    
    ## update z ##
    Z_new = sapply(1:n, function(i){
      #determine whether ith component is a singleton (the only one node in the cluster)
      current.cluster.i = Z_old[i]
      if (counts_Z_old[current.cluster.i] > 1){
        # not a singleton, have |C|+1 choices
        current.counts.noi = counts_Z_old  #current.counts.noi corresponds to |C|
        current.counts.noi[current.cluster.i] = current.counts.noi[current.cluster.i] - 1
        #finding the probs for sampling process
        current.probs = sapply(1:k_old, function(x) {
          Z_star = Z_old
          Z_star[i] = x
          current.prob = log(gamma + current.counts.noi[x]) + loglike(Z_star, Q_old, mu_old, sigma_old, X, i, n)
          return(current.prob)
        })
        Z_star = Z_old
        Z_star[i] = k_old + 1
        current.probs[k_old + 1] = log(gamma) + logmargs(Z_star, X, i, alpha, beta, rou, kappa, delta, xi) + (Vn[k_old + 1] - Vn[k_old])
        current.probs = exp(current.probs - max(current.probs))
        current.probs = current.probs / sum(current.probs)
        
        #choose the cluster number for ith observation
        cluster.i = sample.int(k_old + 1, size = 1, prob = current.probs)
        return(cluster.i)
      } else {
        # a singleton, have |C| choices
        # delete the current cluster
        current.counts.noi = counts_Z_old
        current.counts.noi[current.cluster.i] = current.counts.noi[current.cluster.i] - 1 - gamma
        
        #finding the probs for sampling process
        current.probs = sapply(1:k_old, function(x) {
          Z_star = Z_old
          Z_star[i] = x
          current.prob = log(gamma + current.counts.noi[x]) + loglike(Z_star, Q_old, mu_old, sigma_old, X, i, n)
          return(current.prob)
        })
        Z_star = Z_old
        Z_star[i] = k_old + 1
        current.probs[k_old + 1] = log(gamma) + logmargs(Z_star, X, i, alpha, beta, rou, kappa, delta, xi) + (Vn[k_old] - Vn[k_old - 1])
        current.probs = exp(current.probs - max(current.probs))
        current.probs = current.probs / sum(current.probs)
        
        #choose the cluster number for ith observation
        cluster.i = sample.int(k_old + 1, size = 1, prob = current.probs)
        
        if (cluster.i > k_old) { # if it belongs to a new cluster
          return(current.cluster.i)
        } else { # if it belongs to a previous cluster
          return(cluster.i)
        }
      }
    })
    counts_Z_new = table(as.factor(Z_new))
    k_new = max(Z_new)
    
    if(k_new > k_old){
      Q_star = matrix(0, k_new, k_new)
      Q_star[1:k_old, 1:k_old] = Q
      for (r in (k_old + 1):k_new) {
        for (s in 1:k_new) {
          Q_star[r, s] = rbeta(1, alpha, beta)
          Q_star[s, r] = Q_star[r, s]
        }
      }
      
      sigma_star = matrix(0, k_new, k_new)
      sigma_star[1:k_old, 1:k_old] = sigma
      for (r in (k_old + 1):k_new) {
        for (s in 1:k_new) {
          sigma_star[r, s] = sqrt(rinvgamma(1, delta, xi))
          sigma_star[s, r] = sigma_star[r, s]
        }
      }
      
      mu_star = matrix(0, k_new, k_new)
      mu_star[1:k_old, 1:k_old] = mu
      for (r in (k_old + 1):k_new) {
        for (s in 1:k_new) {
          mu_star[r, s] = rnorm(1, rou, sigma_star[r, s] / sqrt(kappa))
          mu_star[s, r] = mu_star[r, s]
        }
      }
      
      Q = Q_star
      mu = mu_star
      sigma = sigma_star
    }
    
    used_clusters = sort(unique(Z_new))
    Z = as.integer(factor(Z_new, levels = used_clusters))
    counts_Z = table(as.factor(Z))
    k = length(counts_Z)
    
    Q = Q[used_clusters, used_clusters, drop = FALSE]
    mu = mu[used_clusters, used_clusters, drop = FALSE]
    sigma = sigma[used_clusters, used_clusters, drop = FALSE]
    
    ## update A: the adjacency matrix, a n by n matrix ##
    A = matrix(0, n, n)
    for(i in 1:n) {
      r = Z[i]
      for(j in i:n) {
        s = Z[j]
        if(i == j) next
        p_null = (1 - Q[r,s]) * dnorm(X[i,j])
        p_alt = Q[r,s] * dnorm(X[i,j], mu[r,s], sigma[r,s])
        A[i,j] = rbinom(1, size = 1, prob = p_alt / (p_null + p_alt))
      }
    }
    A = A + t(A)
    
    ## update Q, mu, sigma ##
    X_upper = X
    X_upper[lower.tri(X_upper, diag = TRUE)] = NA
    A_upper = A
    A_upper[lower.tri(A_upper, diag = TRUE)] = NA
    
    for (r in 1:k){
      index_r = which(Z == r)
      n_r = length(index_r)
      for (s in r:k)
      {
        index_s = which(Z == s)
        n_s = length(index_s)
        
        alpha.l = 0
        beta.l  = 0
        delta.l = delta
        xi.l    = xi
        rou.l   = rou
        kappa.l = kappa
        
        if (r == s) {
          if (n_r > 1) {
            alpha.l = sum(A_upper[index_r, index_s, drop = FALSE] == 1, na.rm = TRUE)
            beta.l = sum(A_upper[index_r, index_s, drop = FALSE] == 0, na.rm = TRUE)
            
            if(alpha.l > 0) {
              x_values = X_upper[index_r, index_s, drop = FALSE][which(A_upper[index_r, index_s, drop = FALSE] == 1, arr.ind = TRUE)]
              n.l = length(x_values)
              xmean.l = mean(x_values)
              xss.l = sum((x_values - xmean.l) ^ 2)
              
              delta.l = delta + n.l / 2
              xi.l = xi + (xss.l + kappa * n.l / (kappa + n.l) * (xmean.l - rou) ^ 2) / 2
              rou.l = (kappa * rou + n.l * xmean.l) / (kappa + n.l)
              kappa.l = kappa + n.l
            }
          }
        } else {
          if (n_r > 0 && n_s > 0) {
            index_rs = expand.grid(index_r, index_s)
            index_rs = cbind(pmin(index_rs[,1], index_rs[,2]), pmax(index_rs[,1], index_rs[,2]))
            alpha.l = sum(A_upper[index_rs, drop = FALSE] == 1, na.rm = TRUE)
            beta.l = sum(A_upper[index_rs, drop = FALSE] == 0, na.rm = TRUE)
            
            if(alpha.l > 0){
              x_values = X_upper[index_rs, drop = FALSE][which(A_upper[index_rs, drop = FALSE] == 1, arr.ind = TRUE)]
              n.l = length(x_values)
              xmean.l = mean(x_values)
              xss.l = sum((x_values - xmean.l) ^ 2)
              
              delta.l = delta + n.l / 2
              xi.l = xi + (xss.l + kappa * n.l / (kappa + n.l) * (xmean.l - rou) ^ 2) / 2
              rou.l = (kappa * rou + n.l * xmean.l) / (kappa + n.l)
              kappa.l = kappa + n.l
            }
          }
        }
        Q[r,s] = rbeta(1, alpha.l + alpha, beta.l + beta)
        Q[s,r] = Q[r,s]
        sigma[r,s] = sqrt(rinvgamma(1, delta.l, xi.l))
        sigma[s,r] = sigma[r,s]
        mu[r,s] = rnorm(1, rou.l, sigma[r,s] / sqrt(kappa.l))
        mu[s,r] = mu[r,s]
      }
    }
    
    # History[[niter]] = list(Zout = Z, Aout = A)
    if (niter %% 10 == 0) {
      cat("Iteration:", niter, "Cluster Number:", k, "\n", Z,"\n")
    }
  }# for loop over iterations
  
  # calculate l-value
  index = unname(which(upper.tri(X), arr.ind = TRUE)) # the index of all possible edges
  l = matrix(0, n, n)
  l_value = rep(0, n * (n - 1) / 2)
  Z_one_hot = matrix(0, k, n) # the one-hot embedding of cluster assignment
  Z_one_hot[matrix(c(Z, 1:n), ncol = 2)] = 1
  X_value = X[upper.tri(X)]
  density_null = dnorm(X_value, 0, 1) # the probability of data under null hypothesis
  for(r in 1:k){
    for(s in r:k){
      density_alt = dnorm(X_value, mu[r,s], sigma[r,s]) # the probability of data under alternative hypothesis
      if(r == s) {
        mask = (Z_one_hot[r, index[, 1]] * Z_one_hot[s, index[, 2]])
      } else {
        mask = (Z_one_hot[r, index[, 1]] * Z_one_hot[s, index[, 2]]) +
          (Z_one_hot[r, index[, 2]] * Z_one_hot[s, index[, 1]])
      }
      l_value = l_value + density_null * (1 - Q[r,s]) / (density_alt * Q[r,s] + density_null * (1 - Q[r,s])) * mask
    }
  }
  l[upper.tri(l)] = l_value
  l = l + t(l)
  diag(l) = NA
  
  # calculate L-value
  L_function = function(r, s, delta){
    L0 = matrix(0, n, n)
    L1 = matrix(0, n, n)
    L_value = matrix(0, n * (n - 1) / 2, 2)
    L_value[l_value == 1,] = 1
    index_L = (l_value > 0) & (l_value < 1)
    if(length(index_L) > 0){
      a = 1 / sigma[r,s] ^ 2 - 1
      b = -2 * mu[r,s] / sigma[r,s] ^ 2
      c = mu[r,s] ^ 2 / sigma[r,s] ^ 2 + 2 * log(sigma[r,s] * (1 / Q[r,s] - 1) * (1 / l_value - 1))
      if(a != 0){
        L_value[index_L,] = L_value[index_L,] + (a < 0)
        index_L_discriminant = (b ^ 2 > (4 * a * c)) & index_L
        if(sum(index_L_discriminant) > 0){
          d = (-b + matrix(sqrt(pmax(0, b ^ 2 - 4 * a * c[index_L_discriminant])), ncol = 1) %*% c(1, -1)) / (2 * a)
          L_value[index_L_discriminant, 2] = L_value[index_L_discriminant, 2] + pnorm(d[,1], mu[r,s], sigma[r,s]) - pnorm(d[,2], mu[r,s], sigma[r,s])
          L_value[index_L_discriminant, 1] = L_value[index_L_discriminant, 1] + pnorm(d[,1]) - pnorm(d[,2])
        }
      } else {
        if(b != 0) {
          L_value[index_L, 2] = ifelse(b < 0, 
                                       1 - pnorm(-c[index_L] / b, mu[r,s], sigma[r,s]),
                                       pnorm(-c[index_L] / b, mu[r,s], sigma[r,s]))
          L_value[index_L, 1] = ifelse(b < 0, 
                                       1 - pnorm(-c[index_L] / b),
                                       pnorm(-c[index_L] / b))
        } else {
          L_value[index_L,] = 1 * (l_value[index_L] >= 1 - Q[r,s])
        }
      }
    }
    L0[upper.tri(L0)] = L_value[,1]
    L1[upper.tri(L1)] = L_value[,2]
    
    if(delta == 0){
      return(L0)
    } else if (delta == 1){
      return(L1)
    }
  }
  
  # calculate q-value
  q = matrix(0, n, n)
  prob = as.numeric(table(factor(Z, levels = 1:k))) / n
  prob_null = 0
  prob_alt = 0
  for(r in 1:k){
    for(s in r:k){
      q0 = prob[r] * prob[s] * L_function(r, s, delta = 0)
      q1 = prob[r] * prob[s] * L_function(r, s, delta = 1)
      f = 1 + (r != s)
      prob_null = prob_null + f * (1 - Q[r,s]) * q0
      prob_alt = prob_alt + f * Q[r,s] * q1
    }
  }
  index = which(prob_null + prob_alt != 0, arr.ind = TRUE)
  q[index] = prob_null[index] / (prob_alt[index] + prob_null[index])
  q = q + t(q)
  diag(q) = NA
  
  return(list(q = q, Z = Z))
}

MFMSBM_inference = function(q, tau) {
  # calculate phi-value
  n = nrow(q)
  phi = matrix(0, n, n)
  phi[lower.tri(phi)] = (q[lower.tri(q)] < tau)
  phi = phi + t(phi)
  diag(phi) = NA
  return(phi)
}

###### simulation study
## data generation
# NSBM setting
NSBM_generation = function(){
  n = 200
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

## take the data into the MFM-SBM algorithm
MFMSBM_performance = function(data_generation_method, tau = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.25)){
  rindex_MFM = NA
  start = Sys.time()
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
  }
  n = nrow(X)
  end_data = Sys.time()
  elapse_data = end_data - start_data
  cat("Data generation ends: \n")
  print(elapse_data)
  
  cat("MFMSBM fitting begins. \n")
  start_mfm = Sys.time()
  if(data_generation_method == "NSBM"){
    start = Sys.time()  
    niterations = 300
    delta = 2
    xi = 0.5
    rou = 2.5
    kappa = 1  
    alpha = 4
    beta = 7
    gamma = 10
    lambda = 2
    fit_MFM = MFMSBM_estimation(X, niterations, delta, xi, rou, kappa, alpha, beta, gamma, lambda)
    infer_MFM = lapply(tau, function(x) MFMSBM_inference(fit_MFM$q, x))
    TDR_MFM =  sapply(infer_MFM, function(x) sum(A * x, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE))
    FDR_MFM = sapply(infer_MFM, function(x) sum((1 - A) * x, na.rm = TRUE) / pmax(sum(x, na.rm = TRUE), 1))
    rindex_MFM = adj.rand.index(fit_MFM$Z, Z)
  } else if(data_generation_method == "star"){
    niterations = 300
    delta = 2
    xi = 1.5
    rou = 1.5
    kappa = 30
    alpha = 1
    beta = 8
    gamma = 1
    lambda = .4
    fit_MFM = MFMSBM_estimation(X, niterations, delta, xi, rou, kappa, alpha, beta, gamma, lambda)
    infer_MFM = lapply(tau, function(x) MFMSBM_inference(fit_MFM$q, x))
    TDR_MFM =  sapply(infer_MFM, function(x) sum(A * x, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE))
    FDR_MFM = sapply(infer_MFM, function(x) sum((1 - A) * x, na.rm = TRUE) / pmax(sum(x, na.rm = TRUE), 1))
  } else if(data_generation_method == "randombi"){
    niterations = 300
    delta = 2
    xi = 2
    rou = 1
    kappa = 1
    alpha = .5
    beta = 2
    gamma = 3
    lambda = 1
    fit_MFM = MFMSBM_estimation(X, niterations, delta, xi, rou, kappa, alpha, beta, gamma, lambda)
    infer_MFM = lapply(tau, function(x) MFMSBM_inference(fit_MFM$q, x))
    TDR_MFM =  sapply(infer_MFM, function(x) sum(A * x, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE))
    FDR_MFM = sapply(infer_MFM, function(x) sum((1 - A) * x, na.rm = TRUE) / pmax(sum(x, na.rm = TRUE), 1))
  } else if(data_generation_method == "preferattach"){
    niterations = 300
    delta = 5
    xi = .5
    rou = 2
    kappa = 10
    alpha = 1
    beta = 3
    gamma = 2
    lambda = 1
    fit_MFM = MFMSBM_estimation(X, niterations, delta, xi, rou, kappa, alpha, beta, gamma, lambda)
    infer_MFM = lapply(tau, function(x) MFMSBM_inference(fit_MFM$q, x))
    TDR_MFM =  sapply(infer_MFM, function(x) sum(A * x, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE))
    FDR_MFM = sapply(infer_MFM, function(x) sum((1 - A) * x, na.rm = TRUE) / pmax(sum(x, na.rm = TRUE), 1))
  }
  end_mfm = Sys.time()
  elapse_mfm = end_mfm - start_mfm
  cat("MFMSBM fitting ends: \n")
  print(elapse_mfm)
  
  start_rbfk = Sys.time()
  cat("Rebafka's algorithm starts.\n")
  #Rebafka
  fit_rbfk = noisySBM::fitNSBM(X, model = "Gauss01")
  infer_rbfk = lapply(tau, function(x) noisySBM::graphInference(X, nodeClustering = fit_rbfk[[2]]$clustering, theta = fit_rbfk[[2]]$theta, alpha = x, modelFamily = "Gauss")$A)
  TDR_rbfk =  sapply(infer_rbfk, function(x) sum(A * x, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE))
  FDR_rbfk = sapply(infer_rbfk, function(x) sum((1 - A) * x, na.rm = TRUE) / pmax(sum(x, na.rm = TRUE), 1))
  rindex_rbfk = adj.rand.index(fit_rbfk[[2]]$clustering, Z)
  end_rbfk = Sys.time()
  elapse_rbfk = end_rbfk - start_rbfk
  cat("Rebafka's algorithm ends: \n")
  print(elapse_rbfk)
  
  start_klln = Sys.time()
  cat("Kilian's algorithm starts.\n")
  #Kilian
  fit_klln = main_noisySBM_fit(X, NIG = TRUE)
  infer_klln = lapply(tau, function(x) main_noisySBM_infer(fit_klln$Rho, fit_klln$theta, fit_klln$Z, X, alpha = x)$A)
  TDR_klln =  sapply(infer_klln, function(x) sum(A * x, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE))
  FDR_klln = sapply(infer_klln, function(x) sum((1 - A) * x, na.rm = TRUE) / pmax(sum(x, na.rm = TRUE), 1))
  rindex_klln = adj.rand.index(fit_klln$Z, Z)
  end_klln = Sys.time()
  elapse_klln = end_klln - start_klln
  cat("Killian's algorithm ends: \n")
  print(elapse_klln)
  
  end = Sys.time()
  elapse = end - start
  cat("All the models end.\n")
  print(elapse)
  return(list(TDR_MFM = TDR_MFM, FDR_MFM = FDR_MFM, 
              TDR_rbfk = TDR_rbfk, FDR_rbfk = FDR_rbfk, 
              TDR_klln = TDR_klln, FDR_klln = FDR_klln,
              elapse_mfm = elapse_mfm, elapse_rbfk = elapse_rbfk, elapse_klln = elapse_klln,
              rindex_mfm = rindex_MFM, rindex_rbfk = rindex_rbfk, rindex_klln = rindex_klln))
}

cluster_apply = function(iteration, data_generation_method = "NSBM") {
  id = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", "0"))
  set.seed(iteration + 1000 * id)
  result = MFMSBM_performance(data_generation_method)
  return(list(
    TDR_MFM = result$TDR_MFM, FDR_MFM = result$FDR_MFM,
    TDR_rbfk = result$TDR_rbfk, FDR_rbfk = result$FDR_rbfk,
    TDR_klln = result$TDR_klln, FDR_klln = result$FDR_klln,
    elapse_mfm = result$elapse_mfm, 
    elapse_rbfk = result$elapse_rbfk, 
    elapse_klln = result$elapse_klln,
    rindex_mfm = result$rindex_mfm,
    rindex_rbfk = result$rindex_rbfk,
    rindex_klln = result$rindex_klln
  ))
}

MFMSBM_main = function(){
  library(rslurm)
  outdir <- Sys.getenv("OUTDIR", unset = file.path(getwd(), "Output/run_mfmsbm_sim"))
  if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  options(bitmapType = "cairo")
  
  niteration = 500
  cat("Submitting Slurm job...\n")
  sjob = slurm_apply(
    cluster_apply,
    params = data.frame(iteration = 1:niteration),
    jobname = "mfmsbm-sim",
    nodes = niteration,
    cpus_per_node = 1,
    global_objects = c(
      "loglike", "logmargs", "MFMSBM_estimation", "MFMSBM_inference", "MFMSBM_performance",
      "NSBM_generation", "stargraph", "randombipartitegraph", "preferattachgraph", "main_noisySBM_fit", "main_noisySBM_infer"
    ),
    pkgs = c(
      "dplyr","tidyr","ggplot2","fossil","invgamma","MASS",
      "igraph","noisySBM","noisysbmGGM"
    )
  )
  cat("Slurm job submitted!\n")
  results = get_slurm_out(sjob, outtype = "raw")
}

if (sys.nframe() == 0) {
  MFMSBM_main()
}