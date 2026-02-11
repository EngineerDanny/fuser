.l1new_env <- new.env(parent = emptyenv())
.l1new_env$lastNumIters <- NA_integer_
.l1new_env$lastActiveEdges <- NA_integer_
.l1new_env$lastTrace <- NULL
.l1new_env$lastScreening <- NULL

#' Iteration count from the previous fused-L1 new solver run.
#'
#' @return Number of iterations performed in the previous call.
#' @export
fusedLassoProximalNewIterationsTaken <- function() {
  .l1new_env$lastNumIters
}

#' Number of active fusion edges from the previous fused-L1 new solver run.
#'
#' @return Number of active fusion edges in the last fit.
#' @export
fusedLassoProximalNewActiveEdges <- function() {
  .l1new_env$lastActiveEdges
}

#' Iteration diagnostics from the previous fused-L1 new solver run.
#'
#' @return Data frame with per-iteration diagnostics, or NULL if not recorded.
#' @export
fusedLassoProximalNewLastTrace <- function() {
  .l1new_env$lastTrace
}

#' Screening statistics from the previous fused-L1 new solver run.
#'
#' @return A list with screening diagnostics, or NULL if unavailable.
#' @export
fusedLassoProximalNewLastScreening <- function() {
  .l1new_env$lastScreening
}

.clip_unit <- function(x) {
  x[x > 1] <- 1
  x[x < -1] <- -1
  x
}

.validate_l1new_common <- function(lambda, gamma, mu, edge.block) {
  if (!is.numeric(lambda) || length(lambda) != 1L || is.na(lambda) || lambda < 0) {
    stop("Parameter 'lambda' must be a non-negative scalar.")
  }
  if (!is.numeric(gamma) || length(gamma) != 1L || is.na(gamma) || gamma < 0) {
    stop("Parameter 'gamma' must be a non-negative scalar.")
  }
  if (!is.numeric(mu) || length(mu) != 1L || is.na(mu) || mu <= 0) {
    stop("Parameter 'mu' must be a positive scalar.")
  }
  if (!is.numeric(edge.block) || length(edge.block) != 1L || edge.block < 1) {
    stop("Parameter 'edge.block' must be a positive integer.")
  }
}

.init_trace_state <- function(diagnostics = FALSE, trace_every = 1L) {
  e <- new.env(parent = emptyenv())
  e$enabled <- isTRUE(diagnostics)
  e$trace_every <- max(1L, as.integer(trace_every))
  e$iter_global <- 0L
  e$n <- 0L
  e$rows <- list()
  e
}

.append_trace_row <- function(trace_state, row) {
  if (is.null(trace_state) || !isTRUE(trace_state$enabled)) {
    return(invisible(NULL))
  }
  trace_state$n <- trace_state$n + 1L
  trace_state$rows[[trace_state$n]] <- row
  invisible(NULL)
}

.trace_state_to_df <- function(trace_state) {
  if (is.null(trace_state) || !isTRUE(trace_state$enabled) || trace_state$n == 0L) {
    return(NULL)
  }
  rows <- lapply(trace_state$rows, as.data.frame, stringsAsFactors = FALSE)
  do.call(rbind, rows)
}

.objective_l1_fused <- function(B, prep, lambda, gamma, intercept) {
  B.pen <- if (intercept) B[seq_len(prep$p), , drop = FALSE] else B

  sparsity_w <- matrix(prep$penalty.factors, nrow = prep$p, ncol = prep$k)
  sparsity_pen <- lambda * sum(abs(B.pen) * sparsity_w)

  fusion_pen <- 0
  if (gamma > 0 && length(prep$edge.u) > 0) {
    diff.abs <- abs(B.pen[, prep$edge.u, drop = FALSE] - B.pen[, prep$edge.v, drop = FALSE])
    edge_w_mat <- matrix(prep$edge.w, nrow = prep$p, ncol = length(prep$edge.w), byrow = TRUE)
    fusion_pen <- gamma * sum(diff.abs * edge_w_mat)
  }

  loss <- 0
  if (!is.null(prep$XX)) {
    for (g in seq_len(prep$k)) {
      b <- B[, g, drop = FALSE]
      qf <- drop(t(b) %*% prep$XX[[g]] %*% b)
      lin <- drop(t(b) %*% prep$XY[[g]])
      loss <- loss + (qf - 2 * lin + prep$y2[g]) / (2 * prep$samp.sizes[g])
    }
  } else {
    for (g in seq_len(prep$k)) {
      resid <- prep$Y.list[[g]] - prep$X.list[[g]] %*% B[, g, drop = FALSE]
      loss <- loss + sum(resid^2) / (2 * prep$samp.sizes[g])
    }
  }

  loss + sparsity_pen + fusion_pen
}

.compute_delta_lik <- function(XX, XY, X.list, Y, groups, group.names, B.old, samp.sizes) {
  k <- length(group.names)

  if (!is.null(XX)) {
    delta.lik <- lapply(seq_len(k), function(k.i) {
      (XX[[k.i]] %*% B.old[, k.i] - XY[[k.i]]) / samp.sizes[k.i]
    })
  } else {
    delta.lik <- lapply(seq_len(k), function(k.i) {
      g <- group.names[k.i]
      (crossprod(X.list[[k.i]], X.list[[k.i]] %*% B.old[, k.i]) -
         crossprod(X.list[[k.i]], Y[groups == g])) / samp.sizes[k.i]
    })
  }

  do.call("cbind", delta.lik)
}

.init_operator_state <- function(p, k, intercept, B_init = NULL) {
  nrow.beta <- if (intercept) p + 1L else p
  if (is.null(B_init)) {
    B0 <- matrix(0, nrow.beta, k)
  } else {
    B0 <- B_init
  }
  list(
    B.old = B0,
    W = B0,
    weighted.delta.f = matrix(0, nrow(B0), ncol(B0))
  )
}

.extract_B_fusion <- function(B, p, intercept) {
  if (intercept) {
    rbind(B[seq_len(p), , drop = FALSE], 0)
  } else {
    B
  }
}

.apply_group_threshold <- function(B, p, k, intercept, lambda, samp.sizes) {
  B.out <- B
  for (j in seq_len(k)) {
    thresh <- lambda / samp.sizes[j]
    if (intercept) {
      B.out[seq_len(p), j] <- sign(B.out[seq_len(p), j]) *
        pmax(abs(B.out[seq_len(p), j]) - thresh, 0)
    } else {
      B.out[, j] <- sign(B.out[, j]) * pmax(abs(B.out[, j]) - thresh, 0)
    }
  }
  B.out
}

.edge_violation_scores <- function(B, p, intercept, edge.u, edge.v, edge.w, gamma, edge.block = 256L) {
  if (gamma <= 0 || length(edge.u) == 0) {
    return(numeric(0))
  }

  edge.block <- max(1L, as.integer(edge.block))
  B.fusion <- .extract_B_fusion(B, p, intercept)
  num.edges <- length(edge.u)
  scores <- numeric(num.edges)
  block.starts <- seq.int(1L, num.edges, by = edge.block)

  for (start in block.starts) {
    end <- min(start + edge.block - 1L, num.edges)
    idx <- start:end
    diff.block <- B.fusion[, edge.u[idx], drop = FALSE] - B.fusion[, edge.v[idx], drop = FALSE]
    scores[idx] <- apply(abs(diff.block), 2L, max)
  }

  scores * (gamma * edge.w)
}

.max_node_weight <- function(k, edge.u, edge.v, edge.w) {
  if (length(edge.u) == 0) {
    return(0)
  }
  deg.w <- numeric(k)
  for (i in seq_along(edge.u)) {
    deg.w[edge.u[i]] <- deg.w[edge.u[i]] + edge.w[i]
    deg.w[edge.v[i]] <- deg.w[edge.v[i]] + edge.w[i]
  }
  max(deg.w)
}

.fusion_operator_term <- function(
  B.fusion,
  k,
  edge.u,
  edge.v,
  edge.w,
  gamma,
  mu,
  edge.block,
  return_stats = FALSE
) {
  if (gamma <= 0 || length(edge.u) == 0) {
    zero <- matrix(0, nrow(B.fusion), k)
    if (!return_stats) {
      return(zero)
    }
    return(list(
      delta = zero,
      linear = 0,
      prox_sq = 0,
      prox_inf = 0
    ))
  }

  out <- matrix(0, nrow(B.fusion), k)
  num.edges <- length(edge.u)
  edge.block <- max(1L, as.integer(edge.block))
  block.starts <- seq.int(1L, num.edges, by = edge.block)
  linear_sum <- 0
  prox_sq_sum <- 0
  prox_inf <- 0

  for (start in block.starts) {
    end <- min(start + edge.block - 1L, num.edges)
    for (e in start:end) {
      u <- edge.u[e]
      v <- edge.v[e]
      w <- edge.w[e]
      diff <- B.fusion[, u] - B.fusion[, v]
      z <- (gamma * w) * diff
      prox <- pmax(pmin(z / mu, 1), -1)
      contrib <- gamma * w * prox
      out[, u] <- out[, u] + contrib
      out[, v] <- out[, v] - contrib

      if (return_stats) {
        linear_sum <- linear_sum + sum(z * prox)
        prox_sq_sum <- prox_sq_sum + sum(prox^2)
        prox_inf <- max(prox_inf, max(abs(prox)))
      }
    }
  }

  if (!return_stats) {
    return(out)
  }

  list(
    delta = out,
    linear = linear_sum,
    prox_sq = prox_sq_sum,
    prox_inf = prox_inf
  )
}

.run_operator_iterations <- function(
  Y,
  groups,
  group.names,
  XX,
  XY,
  X.list,
  p,
  k,
  samp.sizes,
  lambda,
  gamma,
  mu,
  tol,
  num.it,
  lam.max,
  intercept,
  penalty.factors,
  edge.u,
  edge.v,
  edge.w,
  edge.block,
  state,
  trace_state = NULL,
  stage = "operator",
  prep = NULL
) {
  max.node.weight <- .max_node_weight(k, edge.u, edge.v, edge.w)
  L_U <- lam.max + (lambda^2 + 2 * gamma^2 * max.node.weight) / mu
  L_U.inv <- 1 / L_U

  B.old <- state$B.old
  W <- state$W
  weighted.delta.f <- state$weighted.delta.f

  iter <- 0L
  converged <- FALSE
  for (i in seq_len(num.it)) {
    iter <- i
    if (!is.null(trace_state) && isTRUE(trace_state$enabled)) {
      trace_state$iter_global <- trace_state$iter_global + 1L
    }

    if (!is.null(penalty.factors)) {
      B.sparsity <- B.old[seq_len(p), , drop = FALSE] *
        matrix(penalty.factors, nrow = p, ncol = k)
    } else {
      B.sparsity <- B.old[seq_len(p), , drop = FALSE]
    }

    if (intercept) {
      B.sparsity <- rbind(B.sparsity, 0)
      B.fusion <- rbind(B.old[seq_len(p), , drop = FALSE], 0)
    } else {
      B.fusion <- B.old
    }

    delta.lik <- .compute_delta_lik(
      XX = XX, XY = XY, X.list = X.list, Y = Y, groups = groups,
      group.names = group.names, B.old = B.old, samp.sizes = samp.sizes
    )

    z.sparsity <- lambda * B.sparsity
    A.sparsity <- .clip_unit(z.sparsity / mu)
    sparsity_linear <- sum(z.sparsity * A.sparsity)
    sparsity_sq <- sum(A.sparsity^2)
    sparsity_inf <- max(abs(A.sparsity))

    fusion_stats <- .fusion_operator_term(
      B.fusion = B.fusion, k = k,
      edge.u = edge.u, edge.v = edge.v, edge.w = edge.w,
      gamma = gamma, mu = mu, edge.block = edge.block,
      return_stats = TRUE
    )
    delta.reg <- lambda * A.sparsity + fusion_stats$delta

    delta.f <- delta.lik + delta.reg
    B.new <- W - L_U.inv * delta.f
    weighted.delta.f <- weighted.delta.f + L_U.inv * 0.5 * i * delta.f
    Z <- -weighted.delta.f
    W <- (i * B.new + 2 * Z) / (i + 2)

    improvement <- sum(abs(B.old - B.new))

    if (!is.null(trace_state) && isTRUE(trace_state$enabled)) {
      if (i == 1L || (i %% trace_state$trace_every) == 0L || i == num.it) {
        objective <- NA_real_
        if (!is.null(prep)) {
          objective <- .objective_l1_fused(
            B = B.new,
            prep = prep,
            lambda = lambda,
            gamma = gamma,
            intercept = intercept
          )
        }
        .append_trace_row(trace_state, list(
          stage = stage,
          iter_local = i,
          iter_global = trace_state$iter_global,
          objective = objective,
          dual_linear = sparsity_linear + fusion_stats$linear,
          dual_quadratic = 0.5 * mu * (sparsity_sq + fusion_stats$prox_sq),
          dual_surrogate_reg = (sparsity_linear + fusion_stats$linear) -
            (0.5 * mu * (sparsity_sq + fusion_stats$prox_sq)),
          dual_feas_inf = max(sparsity_inf, fusion_stats$prox_inf),
          update_l1 = improvement,
          kkt_inf = max(abs(delta.f)),
          kkt_fro = sqrt(mean(delta.f^2))
        ))
      }
    }

    if (improvement < tol * p) {
      converged <- TRUE
      break
    }
    B.old <- B.new
  }

  list(
    B = B.new,
    iter = iter,
    converged = converged,
    state = list(B.old = B.new, W = W, weighted.delta.f = weighted.delta.f)
  )
}

.compute_lam_max <- function(XX, X.list) {
  if (!is.null(XX)) {
    return(sum(vapply(XX, function(x) {
      max(eigen(x, symmetric = TRUE, only.values = TRUE)$values)
    }, numeric(1L))))
  }

  sum(vapply(seq_along(X.list), function(g) {
    bigeigen(X.list[[g]])
  }, numeric(1L)))
}

.compute_zero_grad <- function(prep, intercept) {
  k <- prep$k
  p.use <- if (intercept) prep$p else prep$p

  grads <- matrix(0, nrow = p.use, ncol = k)
  if (!is.null(prep$XY)) {
    for (g in seq_len(k)) {
      xy <- as.matrix(prep$XY[[g]])
      grads[, g] <- -xy[seq_len(p.use), 1] / prep$samp.sizes[g]
    }
    return(grads)
  }

  for (g in seq_len(k)) {
    xg <- prep$X.list[[g]][, seq_len(p.use), drop = FALSE]
    yg <- prep$Y.list[[g]]
    grads[, g] <- -as.numeric(crossprod(xg, yg)) / prep$samp.sizes[g]
  }
  grads
}

.screen_edges_grad_zero <- function(
  prep,
  gamma,
  intercept,
  margin = 0,
  max_drop_frac = 1,
  min_keep = 0L
) {
  num.edges <- length(prep$edge.u)
  if (gamma <= 0 || num.edges == 0) {
    keep <- rep(TRUE, num.edges)
    return(list(
      keep = keep,
      scores = numeric(num.edges),
      thresholds = gamma * prep$edge.w,
      screened = 0L
    ))
  }

  margin <- max(0, as.numeric(margin))
  max_drop_frac <- min(1, max(0, as.numeric(max_drop_frac)))
  min_keep <- max(0L, as.integer(min_keep))

  grads <- .compute_zero_grad(prep, intercept = intercept)
  score <- apply(abs(grads[, prep$edge.u, drop = FALSE] - grads[, prep$edge.v, drop = FALSE]), 2L, max)
  threshold <- gamma * prep$edge.w * (1 - margin)

  keep <- !(score <= threshold)

  # Optionally cap how many edges can be dropped by this heuristic.
  max_drop <- floor(max_drop_frac * num.edges)
  drop_idx <- which(!keep)
  if (length(drop_idx) > max_drop) {
    # Keep the highest-score edges among those initially screened.
    to_restore <- drop_idx[order(score[drop_idx], decreasing = TRUE)[seq_len(length(drop_idx) - max_drop)]]
    keep[to_restore] <- TRUE
  }

  # Enforce minimum number of kept edges.
  keep_n <- sum(keep)
  if (keep_n < min_keep && num.edges > 0) {
    need <- min_keep - keep_n
    if (need > 0) {
      cand <- which(!keep)
      restore <- cand[order(score[cand], decreasing = TRUE)[seq_len(min(need, length(cand)))]]
      keep[restore] <- TRUE
    }
  }

  list(
    keep = keep,
    scores = score,
    thresholds = threshold,
    screened = sum(!keep)
  )
}

.prepare_l1new_inputs <- function(
  X,
  Y,
  groups,
  G,
  intercept,
  penalty.factors,
  conserve.memory,
  scaling
) {
  group.names <- sort(unique(groups))
  k <- length(group.names)
  p <- ncol(X)

  if (!is.matrix(G) || nrow(G) != k || ncol(G) != k) {
    stop("G must be a square matrix with dimensions equal to the number of groups.")
  }

  if (is.null(penalty.factors)) {
    penalty.factors <- rep(1, p)
  }

  samp.sizes <- table(groups)[group.names]
  if (!scaling) {
    samp.sizes[] <- 1
  }

  if (intercept) {
    X <- cbind(X, 1)
  }

  if (!conserve.memory) {
    XX.list <- lapply(group.names, function(g) {
      idx <- groups == g
      crossprod(X[idx, , drop = FALSE])
    })
    XY.list <- lapply(group.names, function(g) {
      idx <- groups == g
      crossprod(X[idx, , drop = FALSE], Y[idx])
    })
    X.list <- NULL
  } else {
    XX.list <- NULL
    XY.list <- NULL
    X.list <- lapply(group.names, function(g) {
      X[groups == g, , drop = FALSE]
    })
    X <- NULL
  }

  Y.list <- lapply(group.names, function(g) {
    Y[groups == g]
  })
  y2 <- vapply(Y.list, function(y.g) sum(y.g^2), numeric(1L))

  edges <- which(upper.tri(G) & (G != 0), arr.ind = TRUE)
  edge.u <- if (nrow(edges)) edges[, 1] else integer(0)
  edge.v <- if (nrow(edges)) edges[, 2] else integer(0)
  edge.w <- if (nrow(edges)) G[edges] else numeric(0)

  list(
    Y = Y,
    groups = groups,
    group.names = group.names,
    p = p,
    k = k,
    samp.sizes = as.numeric(samp.sizes),
    penalty.factors = penalty.factors,
    XX = XX.list,
    XY = XY.list,
    X.list = X.list,
    Y.list = Y.list,
    y2 = y2,
    edge.u = edge.u,
    edge.v = edge.v,
    edge.w = edge.w
  )
}

.finalize_l1new_fit <- function(
  fit.raw,
  prep,
  lambda,
  intercept,
  trace_state = NULL,
  screening_info = NULL
) {
  .l1new_env$lastNumIters <- fit.raw$iterations
  .l1new_env$lastActiveEdges <- fit.raw$active_edges
  .l1new_env$lastTrace <- .trace_state_to_df(trace_state)
  .l1new_env$lastScreening <- screening_info
  if (!fit.raw$converged) {
    warning("Reached max iterations without convergence.")
  }

  B.out <- .apply_group_threshold(
    B = fit.raw$B,
    p = prep$p,
    k = prep$k,
    intercept = intercept,
    lambda = lambda,
    samp.sizes = prep$samp.sizes
  )
  colnames(B.out) <- prep$group.names
  B.out
}
