rm(list = ls())
.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")

library(dplyr)
library(tidyr)
library(ggplot2)

SBM_visualization = function(input_dir, output_dir){
  cat("Output job begins...\n")
  files = list.files(input_dir, pattern = "^results_.*\\.RDS$", full.names = TRUE, recursive = TRUE)
  results = do.call(c, lapply(files, readRDS))
  
  tau = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.25)

  TDR_MFM = do.call(rbind, lapply(results, function(x) if (is.list(x)) x[["TDR_MFM"]] else NULL))
  TDR_MFM_mean = apply(TDR_MFM, 2, mean, na.rm = TRUE)
  niteration = nrow(TDR_MFM)
  TDR_MFM_sd = apply(TDR_MFM, 2, sd, na.rm = TRUE)
  TDR_MFM_upper = pmin(TDR_MFM_mean + qnorm(.975) * TDR_MFM_sd / sqrt(niteration), 1)
  TDR_MFM_lower = pmax(TDR_MFM_mean - qnorm(.975) * TDR_MFM_sd / sqrt(niteration), 0)
  
  FDR_MFM = do.call(rbind, lapply(results, function(x) if (is.list(x)) x[["FDR_MFM"]] else NULL))
  FDR_MFM_mean = apply(FDR_MFM, 2, mean, na.rm = TRUE)
  niteration = nrow(FDR_MFM)
  FDR_MFM_sd = apply(FDR_MFM, 2, sd, na.rm = TRUE)
  FDR_MFM_upper = pmin(FDR_MFM_mean + qnorm(.975) * FDR_MFM_sd / sqrt(niteration), 1)
  FDR_MFM_lower = pmax(FDR_MFM_mean - qnorm(.975) * FDR_MFM_sd / sqrt(niteration), 0)
  
  TDR_rbfk = do.call(rbind, lapply(results, function(x) if (is.list(x)) x[["TDR_rbfk"]] else NULL))
  TDR_rbfk_mean = apply(TDR_rbfk, 2, mean, na.rm = TRUE)
  FDR_rbfk = do.call(rbind, lapply(results, function(x) if (is.list(x)) x[["FDR_rbfk"]] else NULL))
  FDR_rbfk_mean = apply(FDR_rbfk, 2, mean, na.rm = TRUE)
  TDR_klln = do.call(rbind, lapply(results, function(x) if (is.list(x)) x[["TDR_klln"]]  else NULL))
  TDR_klln_mean = apply(TDR_klln, 2, mean, na.rm = TRUE)
  FDR_klln = do.call(rbind, lapply(results, function(x) if (is.list(x)) x[["FDR_klln"]]  else NULL))
  FDR_klln_mean = apply(FDR_klln, 2, mean, na.rm = TRUE)
  
  elapse_mfm = sapply(results, function(x) as.numeric(x[["elapse_mfm"]],  units = "mins"))
  elapse_rbfk = sapply(results, function(x) as.numeric(x[["elapse_rbfk"]],  units = "mins"))
  elapse_klln = sapply(results, function(x) as.numeric(x[["elapse_klln"]],  units = "mins"))
  rindex_mfm = sapply(results, function(x) as.numeric(x[["rindex_mfm"]]))
  
  df_results <- data.frame(
    tau   = rep(tau, 2),
    metric= c(rep("TDR", length(tau)), rep("FDR", length(tau))),
    mean  = c(TDR_MFM_mean, FDR_MFM_mean),
    lower = c(TDR_MFM_lower, FDR_MFM_lower),
    upper = c(TDR_MFM_upper, FDR_MFM_upper)
  )
  
  p1 = 
    ggplot(df_results, aes(x = tau, y = mean, color = metric, fill = metric)) +
    geom_line(linewidth = 1.2) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.2, color = NA) +
    geom_point(size = 2) +
    geom_point(aes(x = tau, y = tau), size = 2, shape = 4, color = "red") +
    labs(title = "MFM-SBM Performance Evaluation (NSBM)", x = "tau", y = "Metrics") +
    theme_minimal() +
    scale_color_manual(values = c("TDR" = "#1f77b4", "FDR" = "#ff7f0e")) +
    scale_fill_manual(values = c("TDR" = "#1f77b4", "FDR" = "#ff7f0e")) +
    scale_x_continuous(breaks = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.25)) +
    theme(panel.grid.major = element_blank())
  if(input_dir == "_rslurm_mfmsbmseq"){
    ggsave(file.path(output_dir, "MFM-SBM Evaluation (Sequential).png"), plot = p1, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_mfmsbmsim"){
    ggsave(file.path(output_dir, "MFM-SBM Evaluation (Simultaneous).png"), plot = p1, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_sgdpmmsbmseq"){
    ggsave(file.path(output_dir, "SGDPMM-SBM Evaluation (Sequential).png"), plot = p1, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_sgdpmmsbmsim"){
    ggsave(file.path(output_dir, "SGDPMM-SBM Evaluation (Simultaneous).png"), plot = p1, width = 6, height = 6, dpi = 300)
  }
  
  df_ROC <- data.frame(
    tau       = tau,
    TDR_MFM   = TDR_MFM_mean,
    FDR_MFM   = FDR_MFM_mean,
    TDR_rbfk  = TDR_rbfk_mean,
    FDR_rbfk  = FDR_rbfk_mean,
    TDR_klln  = TDR_klln_mean,
    FDR_klln  = FDR_klln_mean
  ) %>%
    dplyr::select(tau, contains("FDR"), contains("TDR")) %>%
    tidyr::pivot_longer(
      cols = -tau,
      names_to = c(".value", "Method"),
      names_sep = "_"
    ) %>%
    dplyr::mutate(Method = factor(Method, levels = c("MFM", "rbfk", "klln")))
  
  p2 = 
    ggplot(df_ROC, aes(x = FDR, y = TDR, color = Method, shape = Method)) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c("MFM" = "red", "rbfk" = "blue", "klln" = "green"),
                       labels = c("MFM-SBM", "RBFK", "KLLN")) +
    scale_shape_manual(values = c("MFM" = 17, "rbfk" = 16, "klln" = 15),
                       labels = c("MFM-SBM", "RBFK", "KLLN")) +
    labs(title = "ROC Curve (NSBM)", x = "FDR", y = "TDR")
  if(input_dir == "_rslurm_mfmsbmseq"){
    ggsave(file.path(output_dir, "MFMSBM (Sequential) ROC Curve (NSBM).png"), plot = p2, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_mfmsbmsim"){
    ggsave(file.path(output_dir, "MFMSBM (Simultaneous) ROC Curve (NSBM).png"), plot = p2, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_sgdpmmsbmseq"){
    ggsave(file.path(output_dir, "SGDPMMSBM (Sequential) ROC Curve (NSBM).png"), plot = p2, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_sgdpmmsbmsim"){
    ggsave(file.path(output_dir, "SGDPMMSBM (Simultaneous) ROC Curve (NSBM).png"), plot = p2, width = 6, height = 6, dpi = 300)
  }
  
  df_time = rbind(
    data.frame(method = "MFM-SBM", minutes = elapse_mfm),
    data.frame(method = "RBFK",    minutes = elapse_rbfk),
    data.frame(method = "KLLN",    minutes = elapse_klln)
  )
  
  p3 = 
    ggplot(df_time, aes(x = method, y = minutes, fill = method)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, height = 0, size = 0.9, alpha = 0.35) +
    labs(title = "Runtime Comparison", x = NULL, y = "Elapsed time (minutes)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  if(input_dir == "_rslurm_mfmsbmseq"){
    ggsave(file.path(output_dir, "runtime boxplot of MFMSBM (sequential).png"), p3, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_mfmsbmsim"){
    ggsave(file.path(output_dir, "runtime boxplot of MFMSBM (simultaneous).png"), p3, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_sgdpmmsbmseq"){
    ggsave(file.path(output_dir, "runtime boxplot of SGDPMMSBM (sequential).png"), p3, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_sgdpmmsbmsim"){
    ggsave(file.path(output_dir, "runtime boxplot of SGDPMMSBM (simultaneous).png"), p3, width = 6, height = 6, dpi = 300)
  }
  
  df_rand = data.frame(method = "MFM-SBM", RI = rindex_mfm)
  
  p4 = 
    ggplot(df_rand, aes(x = method, y = RI, fill = method)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, height = 0, size = 0.9, alpha = 0.35) +
    coord_cartesian(ylim = c(0, 1)) + 
    labs(title = "Rand Index of MFMSBM (sequential)", x = NULL, y = "Rand Index") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  if(input_dir == "_rslurm_mfmsbmseq"){
    ggsave(file.path(output_dir, "Rand Index of MFMSBM (sequential).png"), p4, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_mfmsbmsim"){
    ggsave(file.path(output_dir, "Rand Index of MFMSBM (simultaneous).png"), p4, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_sgdpmmsbmseq"){
    ggsave(file.path(output_dir, "Rand Index of SGDPMMSBM (sequential).png"), p4, width = 6, height = 6, dpi = 300)
  } else if(input_dir == "_rslurm_sgdpmmsbmsim"){
    ggsave(file.path(output_dir, "Rand Index of SGDPMMSBM (simultaneous).png"), p4, width = 6, height = 6, dpi = 300)
  }
  
  cat("Output job ends!\n")
}
