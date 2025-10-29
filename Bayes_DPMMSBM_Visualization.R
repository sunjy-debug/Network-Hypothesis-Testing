.libPaths("~/R/x86_64-pc-linux-gnu-library/4.4")

library(dplyr)
library(tidyr)
library(ggplot2)
library(viridisLite)

SBM_visualization = function(input_dir, output_dir){
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  cat("Output job begins...\n")
  files = list.files(input_dir, pattern = "^results_.*\\.RDS$", full.names = TRUE, recursive = TRUE)
  results = do.call(c, lapply(files, readRDS))
  
  niteration = 500
  tau = c(0.005, 0.025, 0.05, 0.1, 0.15, 0.25)
  
  TDR_MFMSBM_seq = do.call(rbind, lapply(results, `[[`, "TDR_MFMSBM_seq"))
  FDR_MFMSBM_seq = do.call(rbind, lapply(results, `[[`, "FDR_MFMSBM_seq"))
  TDR_MFMSBM_sim = do.call(rbind, lapply(results, `[[`, "TDR_MFMSBM_sim"))
  FDR_MFMSBM_sim = do.call(rbind, lapply(results, `[[`, "FDR_MFMSBM_sim"))
  TDR_SGDPMMSBM_seq = do.call(rbind, lapply(results, `[[`, "TDR_SGDPMMSBM_seq"))
  FDR_SGDPMMSBM_seq = do.call(rbind, lapply(results, `[[`, "FDR_SGDPMMSBM_seq"))
  TDR_SGDPMMSBM_sim = do.call(rbind, lapply(results, `[[`, "TDR_SGDPMMSBM_sim"))
  FDR_SGDPMMSBM_sim = do.call(rbind, lapply(results, `[[`, "FDR_SGDPMMSBM_sim"))
  
  TDR_MFMSBM_seq_mean = apply(TDR_MFMSBM_seq, 2, mean, na.rm = TRUE)
  TDR_MFMSBM_seq_sd = apply(TDR_MFMSBM_seq, 2, sd, na.rm = TRUE)
  TDR_MFMSBM_seq_upper = pmin(TDR_MFMSBM_seq_mean + qnorm(.975) * TDR_MFMSBM_seq_sd / sqrt(niteration), 1)
  TDR_MFMSBM_seq_lower = pmax(TDR_MFMSBM_seq_mean - qnorm(.975) * TDR_MFMSBM_seq_sd / sqrt(niteration), 0)
  FDR_MFMSBM_seq_mean = apply(FDR_MFMSBM_seq, 2, mean, na.rm = TRUE)
  FDR_MFMSBM_seq_sd = apply(FDR_MFMSBM_seq, 2, sd, na.rm = TRUE)
  FDR_MFMSBM_seq_upper = pmin(FDR_MFMSBM_seq_mean + qnorm(.975) * FDR_MFMSBM_seq_sd / sqrt(niteration), 1)
  FDR_MFMSBM_seq_lower = pmax(FDR_MFMSBM_seq_mean - qnorm(.975) * FDR_MFMSBM_seq_sd / sqrt(niteration), 0)
  TDR_MFMSBM_sim_mean = apply(TDR_MFMSBM_sim, 2, mean, na.rm = TRUE)
  TDR_MFMSBM_sim_sd = apply(TDR_MFMSBM_sim, 2, sd, na.rm = TRUE)
  TDR_MFMSBM_sim_upper = pmin(TDR_MFMSBM_sim_mean + qnorm(.975) * TDR_MFMSBM_sim_sd / sqrt(niteration), 1)
  TDR_MFMSBM_sim_lower = pmax(TDR_MFMSBM_sim_mean - qnorm(.975) * TDR_MFMSBM_sim_sd / sqrt(niteration), 0)
  FDR_MFMSBM_sim_mean = apply(FDR_MFMSBM_sim, 2, mean, na.rm = TRUE)
  FDR_MFMSBM_sim_sd = apply(FDR_MFMSBM_sim, 2, sd, na.rm = TRUE)
  FDR_MFMSBM_sim_upper = pmin(FDR_MFMSBM_sim_mean + qnorm(.975) * FDR_MFMSBM_sim_sd / sqrt(niteration), 1)
  FDR_MFMSBM_sim_lower = pmax(FDR_MFMSBM_sim_mean - qnorm(.975) * FDR_MFMSBM_sim_sd / sqrt(niteration), 0)
  TDR_SGDPMMSBM_seq_mean = apply(TDR_SGDPMMSBM_seq, 2, mean, na.rm = TRUE)
  TDR_SGDPMMSBM_seq_sd = apply(TDR_SGDPMMSBM_seq, 2, sd, na.rm = TRUE)
  TDR_SGDPMMSBM_seq_upper = pmin(TDR_SGDPMMSBM_seq_mean + qnorm(.975) * TDR_SGDPMMSBM_seq_sd / sqrt(niteration), 1)
  TDR_SGDPMMSBM_seq_lower = pmax(TDR_SGDPMMSBM_seq_mean - qnorm(.975) * TDR_SGDPMMSBM_seq_sd / sqrt(niteration), 0)
  FDR_SGDPMMSBM_seq_mean = apply(FDR_SGDPMMSBM_seq, 2, mean, na.rm = TRUE)
  FDR_SGDPMMSBM_seq_sd = apply(FDR_SGDPMMSBM_seq, 2, sd, na.rm = TRUE)
  FDR_SGDPMMSBM_seq_upper = pmin(FDR_SGDPMMSBM_seq_mean + qnorm(.975) * FDR_SGDPMMSBM_seq_sd / sqrt(niteration), 1)
  FDR_SGDPMMSBM_seq_lower = pmax(FDR_SGDPMMSBM_seq_mean - qnorm(.975) * FDR_SGDPMMSBM_seq_sd / sqrt(niteration), 0)
  TDR_SGDPMMSBM_sim_mean = apply(TDR_SGDPMMSBM_sim, 2, mean, na.rm = TRUE)
  TDR_SGDPMMSBM_sim_sd = apply(TDR_SGDPMMSBM_sim, 2, sd, na.rm = TRUE)
  TDR_SGDPMMSBM_sim_upper = pmin(TDR_SGDPMMSBM_sim_mean + qnorm(.975) * TDR_SGDPMMSBM_sim_sd / sqrt(niteration), 1)
  TDR_SGDPMMSBM_sim_lower = pmax(TDR_SGDPMMSBM_sim_mean - qnorm(.975) * TDR_SGDPMMSBM_sim_sd / sqrt(niteration), 0)
  FDR_SGDPMMSBM_sim_mean = apply(FDR_SGDPMMSBM_sim, 2, mean, na.rm = TRUE)
  FDR_SGDPMMSBM_sim_sd = apply(FDR_SGDPMMSBM_sim, 2, sd, na.rm = TRUE)
  FDR_SGDPMMSBM_sim_upper = pmin(FDR_SGDPMMSBM_sim_mean + qnorm(.975) * FDR_SGDPMMSBM_sim_sd / sqrt(niteration), 1)
  FDR_SGDPMMSBM_sim_lower = pmax(FDR_SGDPMMSBM_sim_mean - qnorm(.975) * FDR_SGDPMMSBM_sim_sd / sqrt(niteration), 0)
  
  TDR_rbfk = do.call(rbind, lapply(results, `[[`, "TDR_rbfk"))
  FDR_rbfk = do.call(rbind, lapply(results, `[[`, "FDR_rbfk"))
  TDR_klln = do.call(rbind, lapply(results, `[[`, "TDR_klln"))
  FDR_klln = do.call(rbind, lapply(results, `[[`, "FDR_klln"))
  TDR_rbfk_mean = apply(TDR_rbfk, 2, mean, na.rm = TRUE)
  FDR_rbfk_mean = apply(FDR_rbfk, 2, mean, na.rm = TRUE)
  TDR_klln_mean = apply(TDR_klln, 2, mean, na.rm = TRUE)
  FDR_klln_mean = apply(FDR_klln, 2, mean, na.rm = TRUE)
  
  elapse_mfmsbmseq = sapply(results, function(x) log(as.numeric(x[["elapse_mfmsbmseq"]],  units = "mins")))
  elapse_mfmsbmsim = sapply(results, function(x) log(as.numeric(x[["elapse_mfmsbmsim"]],  units = "mins")))
  elapse_sgdpmmsbmseq = sapply(results, function(x) log(as.numeric(x[["elapse_sgdpmmsbmseq"]],  units = "mins")))
  elapse_sgdpmmsbmsim = sapply(results, function(x) log(as.numeric(x[["elapse_sgdpmmsbmsim"]],  units = "mins")))
  rindex_MFMSBM_seq = sapply(results, function(x) as.numeric(x[["rindex_MFMSBM_seq"]]))
  rindex_MFMSBM_sim = sapply(results, function(x) as.numeric(x[["rindex_MFMSBM_sim"]]))
  rindex_SGDPMMSBM_seq = sapply(results, function(x) as.numeric(x[["rindex_SGDPMMSBM_seq"]]))
  rindex_SGDPMMSBM_sim = sapply(results, function(x) as.numeric(x[["rindex_SGDPMMSBM_sim"]]))
  
  elapse_rbfk = sapply(results, function(x) log(as.numeric(x[["elapse_rbfk"]],  units = "mins")))
  elapse_klln = sapply(results, function(x) log(as.numeric(x[["elapse_klln"]],  units = "mins")))
  rindex_rbfk = sapply(results, function(x) as.numeric(x[["rindex_rbfk"]]))
  rindex_klln = sapply(results, function(x) as.numeric(x[["rindex_klln"]]))
  
  df_ROC <- data.frame(
    tau       = tau,
    TDR_MFMSBM_seq   = TDR_MFMSBM_seq_mean,
    FDR_MFMSBM_seq   = FDR_MFMSBM_seq_mean,
    TDR_MFMSBM_sim   = TDR_MFMSBM_sim_mean,
    FDR_MFMSBM_sim   = FDR_MFMSBM_sim_mean,
    TDR_SGDPMMSBM_seq   = TDR_SGDPMMSBM_seq_mean,
    FDR_SGDPMMSBM_seq   = FDR_SGDPMMSBM_seq_mean,
    TDR_SGDPMMSBM_sim   = TDR_SGDPMMSBM_sim_mean,
    FDR_SGDPMMSBM_sim   = FDR_SGDPMMSBM_sim_mean,
    TDR_rbfk  = TDR_rbfk_mean,
    FDR_rbfk  = FDR_rbfk_mean,
    TDR_klln  = TDR_klln_mean,
    FDR_klln  = FDR_klln_mean
  ) %>%
    dplyr::select(tau, contains("FDR"), contains("TDR")) %>%
    tidyr::pivot_longer(
      cols = -tau,
      names_to = c(".value", "Method"),
      names_pattern = "^(FDR|TDR)_(.+)$"
    ) %>%
    dplyr::mutate(Method = factor(Method, levels = c("MFMSBM_seq", "MFMSBM_sim", "SGDPMMSBM_seq", "SGDPMMSBM_sim", "rbfk", "klln")))
  
  color = viridis(6, option = "viridis", begin = 0.1, end = 0.9)
  p1 = 
    ggplot(df_ROC, aes(x = FDR, y = TDR, color = Method, shape = Method, group = Method)) +
    geom_vline(data = data.frame(x = tau), aes(xintercept = x), linetype = "dashed", color = "grey50", linewidth = 0.6, alpha = 0.7, inherit.aes = FALSE, show.legend = FALSE) +
    geom_line(linewidth = 1) +
    geom_point(size = 3) +
    scale_color_manual(values = c("MFMSBM_seq" = color[1], "MFMSBM_sim" = color[2], "SGDPMMSBM_seq" = color[3], "SGDPMMSBM_sim" = color[4], "rbfk" = color[5], "klln" = color[6]),
                       labels = c("MFMSBM(seq)", "MFMSBM(sim)", "SGDPMMSBM(seq)", "SGDPMMSBM(sim)", "RBFK", "KLLN")) +
    scale_shape_manual(values = c("MFMSBM_seq" = 17, "MFMSBM_sim" = 16, "SGDPMMSBM_seq" = 15, "SGDPMMSBM_sim" = 14, "rbfk" = 13, "klln" = 12),
                       labels = c("MFMSBM(seq)", "MFMSBM(sim)", "SGDPMMSBM(seq)", "SGDPMMSBM(sim)", "RBFK", "KLLN")) +
    labs(title = "ROC Curve", x = "FDR", y = "TDR")
  ggsave(file.path(output_dir, "ROC Curve.png"), plot = p1, width = 6, height = 6, dpi = 300)
  
  df_time = rbind(
    data.frame(method = "MFMSBM(seq)", minutes = elapse_mfmsbmseq),
    data.frame(method = "MFMSBM(sim)", minutes = elapse_mfmsbmsim),
    data.frame(method = "SGDPMMSBM(seq)", minutes = elapse_sgdpmmsbmseq),
    data.frame(method = "SGDPMMSBM(sim)", minutes = elapse_sgdpmmsbmsim),
    data.frame(method = "RBFK", minutes = elapse_rbfk),
    data.frame(method = "KLLN", minutes = elapse_klln)
  )
  
  p2 = 
    ggplot(df_time, aes(x = method, y = minutes, fill = method)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, height = 0, size = 0.9, alpha = 0.35) +
    labs(title = "Runtime Comparison", x = NULL, y = "log elapsed time (minutes)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(file.path(output_dir, "runtime boxplot.png"), plot = p2, width = 6, height = 6, dpi = 300)
  
  df_rand = rbind(
    data.frame(method = "MFMSBM(seq)", RI = rindex_MFMSBM_seq),
    data.frame(method = "MFMSBM(sim)", RI = rindex_MFMSBM_sim),
    data.frame(method = "SGDPMMSBM(seq)", RI = rindex_SGDPMMSBM_seq),
    data.frame(method = "SGDPMMSBM(sim)", RI = rindex_SGDPMMSBM_sim),
    data.frame(method = "RBFK", RI = rindex_rbfk),
    data.frame(method = "KLLN", RI = rindex_klln)
  )
  
  df_prop = df_rand %>%
    dplyr::group_by(method) %>%
    dplyr::summarise(prop_below = mean(RI < 0.75, na.rm = TRUE), .groups = "drop") %>%
    dplyr::mutate(label = sprintf("RI < 0.75: %.1f%%", 100 * prop_below))
  
  p3 = 
    ggplot(df_rand, aes(x = method, y = RI, fill = method)) +
    geom_boxplot(width = 0.55, outlier.shape = NA, alpha = 0.7) +
    geom_jitter(width = 0.15, height = 0, size = 0.9, alpha = 0.35) +
    geom_text(data = df_prop, aes(x = method, y = 0.25, label = label), inherit.aes = FALSE, vjust = 1, size = 3.5) +
    coord_cartesian(ylim = c(0, 1)) + 
    labs(title = "Rand Index", x = NULL, y = "Rand Index") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none")
  ggsave(file.path(output_dir, "Rand Index.png"), plot = p3, width = 6, height = 6, dpi = 300)
  
  cat("Output job ends!\n")
}
