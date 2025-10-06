rm(list=ls())
## We consider undirected graph without selfloops

library(dplyr)
library(tidyr)
library(ggplot2)
library(fossil)
library(invgamma)
library(MASS)
library(igraph)
library(noisySBM)
library(noisysbmGGM)

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
MFMSBM = function(X, niterations, delta, xi, rou, kappa, alpha, beta, gamma, lambda, tau)
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
  
  History = vector("list", niterations)
  
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
    
    History[[niter]] = list(Zout = Z, Aout = A)
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
  
  # calculate phi-value
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
  n = 40
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
  return(list(A = A, X = X))
}

# star A setting
stargraph = function(){
  n = 40
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
  n = 40
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
  n = 40
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
MFM_performance = function(tau, data_generation_method){
  start = Sys.time()
  TDR_MFM = rep(NA, 1)
  FDR_MFM = rep(NA, 1)
  for(itr in 1: 1){
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
      phi_MFM = MFMSBM(X, niterations, delta, xi, rou, kappa, alpha, beta, gamma, lambda, tau)
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
      phi_MFM = MFMSBM(X, niterations, delta, xi, rou, kappa, alpha, beta, gamma, lambda, tau)
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
      phi_MFM = MFMSBM(X, niterations, delta, xi, rou, kappa, alpha, beta, gamma, lambda, tau)
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
      phi_MFM = MFMSBM(X, niterations, delta, xi, rou, kappa, alpha, beta, gamma, lambda, tau)
    }
    TDR_MFM[itr] =  sum(A * phi_MFM, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE)
    FDR_MFM[itr] = sum((1 - A) * phi_MFM, na.rm = TRUE) / pmax(sum(phi_MFM, na.rm = TRUE), 1)
    end_mfm = Sys.time()
    elapse_mfm = end_mfm - start_mfm
    cat("MFMSBM fitting ends: \n")
    print(elapse_mfm)
    
    start_rbfk = Sys.time()
    cat("Rebafka's algorithm starts.\n")
    #Rebafka
    fit_rbfk = noisySBM::fitNSBM(X, model = "Gauss01")
    infer_rbfk = noisySBM::graphInference(X, nodeClustering = fit_rbfk[[2]]$clustering, theta = fit_rbfk[[2]]$theta, alpha = tau, modelFamily = "Gauss")
    phi_rbfk = infer_rbfk$A
    TDR_rbfk =  sum(A * phi_rbfk, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE)
    FDR_rbfk = sum((1 - A) * phi_rbfk, na.rm = TRUE) / pmax(sum(phi_rbfk, na.rm = TRUE), 1)
    end_rbfk = Sys.time()
    elapse_rbfk = end_rbfk - start_rbfk
    cat("Rebafka's algorithm ends: \n")
    print(elapse_rbfk)
    
    start_klln = Sys.time()
    cat("Kilian's algorithm starts.\n")
    #Kilian
    infer_klln = noisysbmGGM::main_noisySBM(X, NIG = TRUE, alpha = tau)
    phi_klln = infer_klln$A
    TDR_klln =  sum(A * phi_klln, na.rm = TRUE) / sum(as.matrix(A), na.rm = TRUE)
    FDR_klln = sum((1 - A) * phi_klln, na.rm = TRUE) / pmax(sum(phi_klln, na.rm = TRUE), 1)
    end_klln = Sys.time()
    elapse_klln = end_klln - start_klln
    cat("Killian's algorithm ends: \n")
    print(elapse_klln)
  }
  end = Sys.time()
  elapse = end - start
  cat("All the models end.\n")
  print(elapse)
  return(list(TDR_MFM = TDR_MFM, FDR_MFM = FDR_MFM, TDR_rbfk = TDR_rbfk, FDR_rbfk = FDR_rbfk, TDR_klln = TDR_klln, FDR_klln = FDR_klln))
}

library(parallel)
cl = makeCluster(detectCores() - 1, type = "PSOCK")
parallel::clusterSetRNGStream(cl, 2025)
tau = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.25)
clusterEvalQ(cl, {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(fossil)
  library(invgamma)
  library(MASS)
  library(igraph)
  library(noisySBM)
  library(noisysbmGGM)
})
clusterExport(cl, c("loglike", "logmargs", "MFMSBM", "MFM_performance", 
                    "NSBM_generation", "stargraph", "randombipartitegraph", "preferattachgraph"))
results = clusterApply(cl, tau, MFM_performance, data_generation_method = "NSBM")
stopCluster(cl)

TDR_MFM = sapply(results, function(x) x$TDR_MFM)
FDR_MFM = sapply(results, function(x) x$FDR_MFM)
if(is.null(dim(TDR_MFM))){
  TDR_MFM_mean = TDR_MFM
  TDR_MFM_sd = rep(0, length(tau))
} else {
  TDR_MFM_mean = apply(TDR_MFM, 2, mean)
  TDR_MFM_sd = apply(TDR_MFM, 2, sd)
}
TDR_MFM_upper = pmin(TDR_MFM_mean + qnorm(.975) * TDR_MFM_sd, 1)
TDR_MFM_lower = pmax(TDR_MFM_mean - qnorm(.975) * TDR_MFM_sd, 0)
if(is.null(dim(FDR_MFM))){
  FDR_MFM_mean = FDR_MFM
  FDR_MFM_sd = rep(0, length(tau))
} else {
  FDR_MFM_mean = apply(FDR_MFM, 2, mean)
  FDR_MFM_sd = apply(FDR_MFM, 2, sd)
}
FDR_MFM_upper = pmin(FDR_MFM_mean + qnorm(.975) * FDR_MFM_sd, 1)
FDR_MFM_lower = pmax(FDR_MFM_mean - qnorm(.975) * FDR_MFM_sd, 0)

TDR_rbfk = sapply(results, function(x) x$TDR_rbfk)
FDR_rbfk = sapply(results, function(x) x$FDR_rbfk)
if(is.null(dim(TDR_rbfk))){
  TDR_rbfk_mean = TDR_rbfk
} else {
  TDR_rbfk_mean = apply(TDR_rbfk, 2, mean)
}
if(is.null(dim(FDR_rbfk))){
  FDR_rbfk_mean = FDR_rbfk
} else {
  FDR_rbfk_mean = apply(FDR_rbfk, 2, mean)
}

TDR_klln = sapply(results, function(x) x$TDR_klln)
FDR_klln = sapply(results, function(x) x$FDR_klln)
if(is.null(dim(TDR_klln))){
  TDR_klln_mean = TDR_klln
} else {
  TDR_klln_mean = apply(TDR_klln, 2, mean)
}
if(is.null(dim(FDR_klln))){
  FDR_klln_mean = FDR_klln
} else {
  FDR_klln_mean = apply(FDR_klln, 2, mean)
}

df_results = data.frame(
  tau = rep(tau, 2),
  metric = c(rep("TDR", length(tau)), rep("FDR", length(tau))),
  mean = c(TDR_MFM_mean, FDR_MFM_mean),
  lower = c(TDR_MFM_lower, FDR_MFM_lower),
  upper = c(TDR_MFM_upper, FDR_MFM_upper)
)

ggplot(df_results, aes(x = tau, y = mean, color = metric, fill = metric)) +
  geom_line(size = 1.2) +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
  geom_point(size = 2) +
  geom_point(aes(x = tau, y = tau), size = 2, shape = 4, color = "red") +
  labs(title = "Simultaneous MFM-SBM Performance Evaluation (NSBM)", x = "tau", y = "Metrics") +
  theme_minimal() +
  scale_color_manual(values = c("TDR" = "#1f77b4", "FDR" = "#ff7f0e")) +
  scale_fill_manual(values = c("TDR" = "#1f77b4", "FDR" = "#ff7f0e")) + 
  scale_x_continuous(breaks = tau) + 
  theme(panel.grid.major = element_blank())
ggsave("Simultaneous MFM-SBM Performance Evaluation (NSBM).png", width = 4, height = 6)

df_ROC = data.frame(
  TDR_MFM = TDR_MFM_mean,
  FDR_MFM = FDR_MFM_mean,
  TDR_rbfk = TDR_rbfk_mean,
  FDR_rbfk = FDR_rbfk_mean,
  TDR_klln = TDR_klln_mean,
  FDR_klln = FDR_klln_mean,
  tau = tau)%>%
  dplyr::select(tau, contains("FDR"), contains("TDR")) %>%
  tidyr::pivot_longer(
    cols = -tau,
    names_to = c(".value", "Method"),
    names_sep = "_") %>%
  dplyr::mutate(Method = factor(Method, levels = c("MFM", "rbfk", "klln")))

ggplot(df_ROC, aes(x = FDR, y = TDR, color = Method, shape = Method)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  scale_color_manual(
    values = c("MFM" = "red", "rbfk" = "blue", "klln" = "green"),
    labels = c("MFM-SBM", "RBFK", "KLLN")
  ) +
  scale_shape_manual(
    values = c("MFM" = 17, "rbfk" = 16, "klln" = 15),
    labels = c("MFM-SBM", "RBFK", "KLLN")
  ) +
  labs(title = "ROC Curve with Simultaneous MFMSBM (NSBM)", x = "FDR", y = "TDR")
ggsave("ROC Curve with Simultaneous MFMSBM (NSBM).png", width = 4, height = 6)
