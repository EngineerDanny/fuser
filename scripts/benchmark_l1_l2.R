#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(fuser))

# -----------------------------
# User defaults (edit here)
# -----------------------------
# Run with:
#   Rscript scripts/benchmark_l1_l2.R
# Optional CLI flags can still override these defaults.
DEFAULTS <- list(
  mode = "baseline",                  # baseline | sweep | run
  method = NULL,                       # used only when mode == "run"
  output_csv = "benchmark_l1_l2_summary.csv",
  seed = 20260206L,
  k = 6L,
  p = 800L,
  n_group_train = 100L,
  n_group_test = 50L,
  sigma = 0.05,
  lambda = 1e-3,
  gamma = 1e-3,
  scaling = FALSE,
  tol = 9e-5,
  num_it = 2000L,
  c_flag = TRUE,
  reps = 3L,
  lambdas = "1e-4,3e-4,1e-3,3e-3,1e-2,3e-2,1e-1",
  gammas = "1e-4,1e-3,1e-2,1e-1"
)

parse_args <- function(argv) {
  out <- list()
  i <- 1
  while (i <= length(argv)) {
    key <- argv[[i]]
    if (!startsWith(key, "--")) {
      stop(sprintf("Unexpected argument: %s", key))
    }
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
  if (is.null(x)) {
    return(default)
  }
  tolower(x) %in% c("true", "t", "1", "yes", "y")
}

as_num_vec <- function(x) {
  if (is.null(x) || nchar(x) == 0) {
    return(numeric(0))
  }
  parts <- trimws(strsplit(x, ",", fixed = TRUE)[[1]])
  parts <- parts[nzchar(parts)]
  vals <- suppressWarnings(as.numeric(parts))
  if (anyNA(vals)) {
    stop(sprintf("Could not parse numeric vector from: %s", x))
  }
  vals
}

predict_from_beta <- function(beta, X, groups) {
  group_levels <- sort(unique(groups))
  if (is.null(colnames(beta))) {
    colnames(beta) <- as.character(group_levels)
  }
  yhat <- numeric(nrow(X))
  for (g in group_levels) {
    idx <- which(groups == g)
    g_col <- match(as.character(g), colnames(beta))
    if (is.na(g_col)) {
      g_col <- which(group_levels == g)
    }
    yhat[idx] <- as.numeric(X[idx, , drop = FALSE] %*% beta[, g_col, drop = FALSE])
  }
  yhat
}

mse <- function(y, yhat) mean((y - yhat)^2)
rmse <- function(y, yhat) sqrt(mse(y, yhat))
mae <- function(y, yhat) mean(abs(y - yhat))
r2 <- function(y, yhat) 1 - sum((y - yhat)^2) / sum((y - mean(y))^2)

generate_data <- function(seed, k, p, n_group_train, n_group_test, sigma) {
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

  y_train <- unlist(lapply(seq_len(k), function(i) {
    as.numeric(X_train_list[[i]] %*% beta_true[, i] + rnorm(n_group_train, 0, sigma))
  }), use.names = FALSE)

  y_test <- unlist(lapply(seq_len(k), function(i) {
    as.numeric(X_test_list[[i]] %*% beta_true[, i] + rnorm(n_group_test, 0, sigma))
  }), use.names = FALSE)

  list(
    X_train = do.call(rbind, X_train_list),
    y_train = y_train,
    groups_train = groups_train,
    X_test = do.call(rbind, X_test_list),
    y_test = y_test,
    groups_test = groups_test,
    beta_true = beta_true,
    G = matrix(1, k, k),
    meta = list(
      seed = seed,
      k = k,
      p = p,
      n_train = length(groups_train),
      n_test = length(groups_test),
      sigma = sigma
    )
  )
}

fit_once <- function(method, data, lambda, gamma, scaling, tol, num_it, c_flag) {
  X_train <- data$X_train
  y_train <- data$y_train
  groups_train <- data$groups_train
  X_test <- data$X_test
  y_test <- data$y_test
  groups_test <- data$groups_test
  beta_true <- data$beta_true
  G <- data$G

  warnings_seen <- character(0)
  fit <- NULL

  timed <- system.time({
    fit <- withCallingHandlers({
      if (method == "l1") {
        fusedLassoProximal(
          X_train, y_train, groups_train,
          lambda = lambda,
          gamma = gamma,
          G = G,
          tol = tol,
          num.it = num_it,
          intercept = FALSE,
          scaling = scaling,
          c.flag = c_flag
        )
      } else if (method == "l2") {
        fusedL2DescentGLMNet(
          X_train, y_train, groups_train,
          lambda = lambda,
          G = G,
          gamma = gamma,
          scaling = scaling
        )
      } else {
        stop(sprintf("Unknown method: %s", method))
      }
    }, warning = function(w) {
      warnings_seen <<- c(warnings_seen, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  })

  yhat_train <- predict_from_beta(fit, X_train, groups_train)
  yhat_test <- predict_from_beta(fit, X_test, groups_test)

  list(
    method = method,
    elapsed = unname(timed[["elapsed"]]),
    user = unname(timed[["user.self"]]),
    sys = unname(timed[["sys.self"]]),
    train_rmse = rmse(y_train, yhat_train),
    train_mae = mae(y_train, yhat_train),
    train_r2 = r2(y_train, yhat_train),
    test_rmse = rmse(y_test, yhat_test),
    test_mae = mae(y_test, yhat_test),
    test_r2 = r2(y_test, yhat_test),
    beta_rmse = sqrt(mean((fit - beta_true)^2)),
    warnings = paste(unique(warnings_seen), collapse = " | ")
  )
}

print_result_line <- function(x, rep_id = NA_integer_) {
  cat(sprintf(
    paste0(
      "RESULT method=%s rep=%s elapsed=%.6f user=%.6f sys=%.6f ",
      "train_rmse=%.6f train_mae=%.6f train_r2=%.6f ",
      "test_rmse=%.6f test_mae=%.6f test_r2=%.6f beta_rmse=%.6f warnings=%s\n"
    ),
    x$method,
    ifelse(is.na(rep_id), "NA", as.character(rep_id)),
    x$elapsed,
    x$user,
    x$sys,
    x$train_rmse,
    x$train_mae,
    x$train_r2,
    x$test_rmse,
    x$test_mae,
    x$test_r2,
    x$beta_rmse,
    ifelse(nchar(x$warnings) == 0, "none", shQuote(x$warnings))
  ))
}

argv <- parse_args(commandArgs(trailingOnly = TRUE))

get_arg <- function(x, key) {
  if (key %in% names(x)) x[[key]] else NULL
}

get_opt <- function(key) {
  cli <- get_arg(argv, key)
  if (!is.null(cli)) {
    return(cli)
  }
  DEFAULTS[[key]]
}

mode <- get_opt("mode")
method <- get_opt("method")
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
tol <- as.numeric(get_opt("tol"))
num_it <- as.integer(get_opt("num_it"))
c_flag <- as_bool(get_opt("c_flag"), default = TRUE)
reps <- as.integer(get_opt("reps"))

lambdas <- as_num_vec(get_opt("lambdas"))
gammas <- as_num_vec(get_opt("gammas"))

d <- generate_data(seed, k, p, n_group_train, n_group_test, sigma)

if (mode == "run") {
  if (is.null(method)) {
    stop("--method is required in mode=run")
  }
  res <- fit_once(method, d, lambda, gamma, scaling, tol, num_it, c_flag)
  print_result_line(res)

  run_df <- data.frame(
    row_type = "run",
    method = res$method,
    rep = 1L,
    gamma = gamma,
    lambda = lambda,
    elapsed = res$elapsed,
    user = res$user,
    sys = res$sys,
    train_rmse = res$train_rmse,
    train_mae = res$train_mae,
    train_r2 = res$train_r2,
    test_rmse = res$test_rmse,
    test_mae = res$test_mae,
    test_r2 = res$test_r2,
    beta_rmse = res$beta_rmse,
    warnings = ifelse(nchar(res$warnings) == 0, "none", res$warnings),
    stringsAsFactors = FALSE
  )
  write.csv(run_df, output_csv, row.names = FALSE)
  cat(sprintf("Saved summary CSV: %s\n", output_csv))
  quit(save = "no", status = 0)
}

if (mode == "baseline") {
  rows <- list()
  idx <- 1L
  for (m in c("l1", "l2")) {
    for (r in seq_len(reps)) {
      res <- fit_once(m, d, lambda, gamma, scaling, tol, num_it, c_flag)
      print_result_line(res, rep_id = r)
      rows[[idx]] <- data.frame(
        row_type = "run",
        method = res$method,
        rep = r,
        gamma = gamma,
        lambda = lambda,
        elapsed = res$elapsed,
        user = res$user,
        sys = res$sys,
        train_rmse = res$train_rmse,
        train_mae = res$train_mae,
        train_r2 = res$train_r2,
        test_rmse = res$test_rmse,
        test_mae = res$test_mae,
        test_r2 = res$test_r2,
        beta_rmse = res$beta_rmse,
        warnings = ifelse(nchar(res$warnings) == 0, "none", res$warnings),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1L
    }
  }

  run_df <- do.call(rbind, rows)
  mean_df <- aggregate(cbind(elapsed, user, sys, train_rmse, train_mae, train_r2, test_rmse, test_mae, test_r2, beta_rmse) ~ method, run_df, mean)
  mean_df$row_type <- "mean"
  mean_df$rep <- NA_integer_
  mean_df$gamma <- gamma
  mean_df$lambda <- lambda
  mean_df$warnings <- "n/a"

  out_df <- rbind(
    run_df[, c("row_type", "method", "rep", "gamma", "lambda", "elapsed", "user", "sys", "train_rmse", "train_mae", "train_r2", "test_rmse", "test_mae", "test_r2", "beta_rmse", "warnings")],
    mean_df[, c("row_type", "method", "rep", "gamma", "lambda", "elapsed", "user", "sys", "train_rmse", "train_mae", "train_r2", "test_rmse", "test_mae", "test_r2", "beta_rmse", "warnings")]
  )

  write.csv(out_df, output_csv, row.names = FALSE)
  cat("\nMean metrics by method:\n")
  print(mean_df[, c("method", "elapsed", "train_rmse", "test_rmse", "test_r2", "beta_rmse")], row.names = FALSE)
  cat(sprintf("\nSaved summary CSV: %s\n", output_csv))
  quit(save = "no", status = 0)
}

if (mode == "sweep") {
  rows <- list()
  idx <- 1L
  for (m in c("l1", "l2")) {
    for (ga in gammas) {
      for (la in lambdas) {
        res <- fit_once(m, d, la, ga, scaling, tol, num_it, c_flag)
        rows[[idx]] <- data.frame(
          row_type = "sweep",
          method = m,
          rep = NA_integer_,
          gamma = ga,
          lambda = la,
          elapsed = res$elapsed,
          user = res$user,
          sys = res$sys,
          train_rmse = res$train_rmse,
          train_mae = res$train_mae,
          train_r2 = res$train_r2,
          test_rmse = res$test_rmse,
          test_mae = res$test_mae,
          test_r2 = res$test_r2,
          beta_rmse = res$beta_rmse,
          warnings = ifelse(nchar(res$warnings) == 0, "none", res$warnings),
          stringsAsFactors = FALSE
        )
        idx <- idx + 1L
        cat(sprintf("SWEEP method=%s gamma=%g lambda=%g test_rmse=%.6f test_r2=%.6f\n",
                    m, ga, la, res$test_rmse, res$test_r2))
      }
    }
  }

  sweep_df <- do.call(rbind, rows)
  best_rows <- do.call(rbind, lapply(split(sweep_df, sweep_df$method), function(sub) {
    sub[which.min(sub$test_rmse), , drop = FALSE]
  }))
  best_rows$row_type <- "best"

  out_df <- rbind(sweep_df, best_rows)
  write.csv(out_df, output_csv, row.names = FALSE)

  cat("\nBest by test RMSE (per method):\n")
  print(best_rows[, c("method", "gamma", "lambda", "elapsed", "test_rmse", "test_r2", "beta_rmse", "warnings")], row.names = FALSE)
  cat(sprintf("\nSaved summary CSV: %s\n", output_csv))
  quit(save = "no", status = 0)
}

stop(sprintf("Unknown mode: %s", mode))
