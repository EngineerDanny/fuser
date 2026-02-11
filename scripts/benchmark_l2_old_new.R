#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(fuser))
source("R/l2_fusion_new.R")

# -----------------------------
# User defaults (edit here)
# -----------------------------
DEFAULTS <- list(
  mode = "baseline",                     # baseline | run
  variant = NULL,                        # used only when mode == "run": old | new
  output_csv = "benchmark_l2_old_new_summary.csv",
  seed = 20260206L,
  k = 40L,
  p = 100L,
  n_group_train = 50L,
  n_group_test = 25L,
  sigma = 0.05,
  lambda = 1e-3,
  gamma = 1e-3,
  scaling = FALSE,
  reps = 1L,
  g_structure = "sparse_chain"           # sparse_chain | dense
)

parse_args <- function(argv) {
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) stop(sprintf("Unexpected argument: %s", key))
    key <- substring(key, 3)
    if (i == length(argv) || startsWith(argv[[i + 1]], "--")) {
      out[[key]] <- "TRUE"
      i <- i + 1
    } else {
      out[[key]] <- argv[[i + 1]]
      i <- i + 2
    }
  }
  out
}

as_bool <- function(x, default = FALSE) {
  if (is.null(x)) return(default)
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

predict_from_beta <- function(beta, X, groups) {
  group_levels <- sort(unique(groups))
  if (is.null(colnames(beta))) colnames(beta) <- as.character(group_levels)
  yhat <- numeric(nrow(X))
  for (g in group_levels) {
    idx <- which(groups == g)
    g_col <- match(as.character(g), colnames(beta))
    if (is.na(g_col)) g_col <- which(group_levels == g)
    yhat[idx] <- as.numeric(X[idx, , drop = FALSE] %*% beta[, g_col, drop = FALSE])
  }
  yhat
}

rmse <- function(y, yhat) sqrt(mean((y - yhat)^2))
r2 <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

build_G <- function(k, structure) {
  if (structure == "dense") {
    return(matrix(1, k, k))
  }
  if (structure == "sparse_chain") {
    G <- matrix(0, k, k)
    if (k > 1) {
      for (i in 1:(k - 1)) {
        G[i, i + 1] <- 1
        G[i + 1, i] <- 1
      }
    }
    return(G)
  }
  stop("Unknown g_structure: ", structure)
}

generate_data <- function(seed, k, p, n_group_train, n_group_test, sigma, g_structure) {
  set.seed(seed)
  groups_train <- rep(seq_len(k), each = n_group_train)
  groups_test <- rep(seq_len(k), each = n_group_test)

  beta_true <- matrix(0, nrow = p, ncol = k)
  nonzero_ind <- rbinom(p * k, 1, 0.02 / k)
  nonzero_shared <- rbinom(p, 1, 0.02)
  beta_true[which(nonzero_ind == 1)] <- rnorm(sum(nonzero_ind), 1, 0.25)
  beta_true[which(nonzero_shared == 1), ] <- rnorm(sum(nonzero_shared), -1, 0.25)

  X_train_list <- lapply(seq_len(k), function(i) matrix(rnorm(n_group_train * p), n_group_train, p))
  X_test_list <- lapply(seq_len(k), function(i) matrix(rnorm(n_group_test * p), n_group_test, p))
  X_train <- do.call(rbind, X_train_list)
  X_test <- do.call(rbind, X_test_list)

  y_train <- unlist(lapply(seq_len(k), function(i) {
    as.numeric(X_train_list[[i]] %*% beta_true[, i] + rnorm(n_group_train, 0, sigma))
  }), use.names = FALSE)
  y_test <- unlist(lapply(seq_len(k), function(i) {
    as.numeric(X_test_list[[i]] %*% beta_true[, i] + rnorm(n_group_test, 0, sigma))
  }), use.names = FALSE)

  list(
    X_train = X_train,
    y_train = y_train,
    groups_train = groups_train,
    X_test = X_test,
    y_test = y_test,
    groups_test = groups_test,
    beta_true = beta_true,
    G = build_G(k, g_structure)
  )
}

fit_once <- function(variant, data, lambda, gamma, scaling) {
  warnings_seen <- character(0)
  fit <- NULL

  t <- system.time({
    fit <- withCallingHandlers({
      if (variant == "old") {
        fuser::fusedL2DescentGLMNet(
          data$X_train, data$y_train, data$groups_train,
          lambda = lambda, G = data$G, gamma = gamma, scaling = scaling
        )
      } else if (variant == "new") {
        fusedL2DescentGLMNetNew(
          data$X_train, data$y_train, data$groups_train,
          lambda = lambda, G = data$G, gamma = gamma, scaling = scaling
        )
      } else {
        stop("Unknown variant: ", variant)
      }
    }, warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  })

  yhat_train <- predict_from_beta(fit, data$X_train, data$groups_train)
  yhat_test <- predict_from_beta(fit, data$X_test, data$groups_test)

  list(
    variant = variant,
    elapsed = unname(t[["elapsed"]]),
    user = unname(t[["user.self"]]),
    sys = unname(t[["sys.self"]]),
    train_rmse = rmse(data$y_train, yhat_train),
    test_rmse = rmse(data$y_test, yhat_test),
    test_r2 = r2(data$y_test, yhat_test),
    beta_rmse = sqrt(mean((fit - data$beta_true)^2)),
    warnings = ifelse(length(warnings_seen) == 0, "none", paste(unique(warnings_seen), collapse = " | "))
  )
}

print_line <- function(x, rep_id, k, p, g_structure) {
  cat(sprintf(
    paste0(
      "RESULT variant=%s rep=%d k=%d p=%d g_structure=%s ",
      "elapsed=%.6f user=%.6f sys=%.6f train_rmse=%.6f test_rmse=%.6f ",
      "test_r2=%.6f beta_rmse=%.6f warnings=%s\n"
    ),
    x$variant, rep_id, k, p, g_structure,
    x$elapsed, x$user, x$sys, x$train_rmse, x$test_rmse,
    x$test_r2, x$beta_rmse, shQuote(x$warnings)
  ))
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))
get_arg <- function(x, key) if (key %in% names(x)) x[[key]] else NULL
get_opt <- function(key) {
  cli <- get_arg(argv, key)
  if (!is.null(cli)) return(cli)
  DEFAULTS[[key]]
}

mode <- get_opt("mode")
variant <- get_opt("variant")
output_csv <- get_opt("output_csv")
seed <- as.integer(get_opt("seed"))
k <- as.integer(get_opt("k"))
p <- as.integer(get_opt("p"))
n_group_train <- as.integer(get_opt("n_group_train"))
n_group_test <- as.integer(get_opt("n_group_test"))
sigma <- as.numeric(get_opt("sigma"))
lambda <- as.numeric(get_opt("lambda"))
gamma <- as.numeric(get_opt("gamma"))
scaling <- as_bool(get_opt("scaling"), default = FALSE)
reps <- as.integer(get_opt("reps"))
g_structure <- as.character(get_opt("g_structure"))

data <- generate_data(seed, k, p, n_group_train, n_group_test, sigma, g_structure)

if (mode == "run") {
  if (is.null(variant)) stop("--variant is required for mode=run")
  res <- fit_once(variant, data, lambda, gamma, scaling)
  print_line(res, 1L, k, p, g_structure)
  run_df <- data.frame(
    row_type = "run",
    variant = res$variant,
    rep = 1L,
    k = k,
    p = p,
    g_structure = g_structure,
    elapsed = res$elapsed,
    user = res$user,
    sys = res$sys,
    train_rmse = res$train_rmse,
    test_rmse = res$test_rmse,
    test_r2 = res$test_r2,
    beta_rmse = res$beta_rmse,
    warnings = res$warnings,
    stringsAsFactors = FALSE
  )
  write.csv(run_df, output_csv, row.names = FALSE)
  cat(sprintf("Saved summary CSV: %s\n", output_csv))
  quit(save = "no", status = 0)
}

if (mode == "baseline") {
  rows <- list()
  idx <- 1L
  for (v in c("old", "new")) {
    for (r in seq_len(reps)) {
      res <- fit_once(v, data, lambda, gamma, scaling)
      print_line(res, r, k, p, g_structure)
      rows[[idx]] <- data.frame(
        row_type = "run",
        variant = res$variant,
        rep = r,
        k = k,
        p = p,
        g_structure = g_structure,
        elapsed = res$elapsed,
        user = res$user,
        sys = res$sys,
        train_rmse = res$train_rmse,
        test_rmse = res$test_rmse,
        test_r2 = res$test_r2,
        beta_rmse = res$beta_rmse,
        warnings = res$warnings,
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }
  run_df <- do.call(rbind, rows)
  mean_df <- aggregate(cbind(elapsed, user, sys, train_rmse, test_rmse, test_r2, beta_rmse) ~ variant, run_df, mean)
  mean_df$row_type <- "mean"
  mean_df$rep <- NA_integer_
  mean_df$k <- k
  mean_df$p <- p
  mean_df$g_structure <- g_structure
  mean_df$warnings <- "n/a"
  out_df <- rbind(
    run_df[, c("row_type", "variant", "rep", "k", "p", "g_structure", "elapsed", "user", "sys", "train_rmse", "test_rmse", "test_r2", "beta_rmse", "warnings")],
    mean_df[, c("row_type", "variant", "rep", "k", "p", "g_structure", "elapsed", "user", "sys", "train_rmse", "test_rmse", "test_r2", "beta_rmse", "warnings")]
  )
  write.csv(out_df, output_csv, row.names = FALSE)
  cat("\nMean metrics by variant:\n")
  print(mean_df[, c("variant", "elapsed", "train_rmse", "test_rmse", "test_r2", "beta_rmse")], row.names = FALSE)
  cat(sprintf("\nSaved summary CSV: %s\n", output_csv))
  quit(save = "no", status = 0)
}

stop("Unknown mode: ", mode)
