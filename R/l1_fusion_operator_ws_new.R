.genFusedLassoProximalOperatorWS <- function(
  prep,
  lambda,
  gamma,
  mu,
  tol,
  num.it,
  lam.max,
  intercept,
  edge.block,
  ws_init_edges,
  ws_add_edges,
  ws_max_outer,
  ws_inner_it,
  ws_violation_tol,
  ws_final_full,
  ws_final_it,
  screening,
  screen_margin,
  screen_max_drop_frac,
  screen_min_keep,
  trace_state = NULL
) {
  if (is.null(lam.max)) {
    lam.max <- .compute_lam_max(prep$XX, prep$X.list)
  }

  if (gamma <= 0 || length(prep$edge.u) == 0) {
    return(.genFusedLassoProximalOperatorOnly(
      prep = prep,
      lambda = lambda,
      gamma = gamma,
      mu = mu,
      tol = tol,
      num.it = num.it,
      lam.max = lam.max,
      intercept = intercept,
      edge.block = edge.block,
      trace_state = trace_state
    ))
  }

  screen_info <- list(
    method = screening,
    total_edges = length(prep$edge.u),
    kept_edges = length(prep$edge.u),
    screened_edges = 0L,
    kept_fraction = if (length(prep$edge.u) > 0) 1 else NA_real_,
    margin = screen_margin,
    max_drop_frac = screen_max_drop_frac,
    min_keep = as.integer(screen_min_keep)
  )

  edge.u <- prep$edge.u
  edge.v <- prep$edge.v
  edge.w <- prep$edge.w
  if (screening == "grad_zero" && gamma > 0 && length(edge.u) > 0) {
    screened <- .screen_edges_grad_zero(
      prep = prep,
      gamma = gamma,
      intercept = intercept,
      margin = screen_margin,
      max_drop_frac = screen_max_drop_frac,
      min_keep = screen_min_keep
    )
    edge.u <- edge.u[screened$keep]
    edge.v <- edge.v[screened$keep]
    edge.w <- edge.w[screened$keep]
    screen_info$kept_edges <- length(edge.u)
    screen_info$screened_edges <- screened$screened
    screen_info$kept_fraction <- if (screen_info$total_edges > 0) {
      screen_info$kept_edges / screen_info$total_edges
    } else {
      NA_real_
    }
  }

  num.edges <- length(edge.u)
  ws_init_edges <- min(as.integer(ws_init_edges), num.edges)
  ws_add_edges <- max(1L, as.integer(ws_add_edges))
  ws_max_outer <- max(1L, as.integer(ws_max_outer))
  if (is.null(ws_inner_it)) {
    ws_inner_it <- max(30L, ceiling(num.it / ws_max_outer))
  }
  ws_inner_it <- max(1L, as.integer(ws_inner_it))
  if (is.null(ws_final_it)) {
    ws_final_it <- num.it
  }
  ws_final_it <- max(1L, as.integer(ws_final_it))

  active <- rep(FALSE, num.edges)
  if (ws_init_edges > 0) {
    init.order <- order(edge.w, decreasing = TRUE)
    active[init.order[seq_len(ws_init_edges)]] <- TRUE
  }

  state <- .init_operator_state(prep$p, prep$k, intercept)
  total.iter <- 0L
  converged <- FALSE

  for (outer in seq_len(ws_max_outer)) {
    idx <- which(active)
    core <- .run_operator_iterations(
      Y = prep$Y,
      groups = prep$groups,
      group.names = prep$group.names,
      XX = prep$XX,
      XY = prep$XY,
      X.list = prep$X.list,
      p = prep$p,
      k = prep$k,
      samp.sizes = prep$samp.sizes,
      lambda = lambda,
      gamma = gamma,
      mu = mu,
      tol = tol,
      num.it = ws_inner_it,
      lam.max = lam.max,
      intercept = intercept,
      penalty.factors = prep$penalty.factors,
      edge.u = edge.u[idx],
      edge.v = edge.v[idx],
      edge.w = edge.w[idx],
      edge.block = edge.block,
      state = state,
      trace_state = trace_state,
      stage = sprintf("ws_outer_%d", outer),
      prep = prep
    )
    total.iter <- total.iter + core$iter
    state <- core$state
    B.cur <- core$B

    scores <- .edge_violation_scores(
      B = B.cur,
      p = prep$p,
      intercept = intercept,
      edge.u = edge.u,
      edge.v = edge.v,
      edge.w = edge.w,
      gamma = gamma,
      edge.block = edge.block
    )
    inactive <- which(!active)
    violating <- inactive[scores[inactive] > ws_violation_tol]
    if (length(violating) == 0) {
      converged <- core$converged
      break
    }

    add.n <- min(ws_add_edges, length(violating))
    add.idx <- violating[order(scores[violating], decreasing = TRUE)[seq_len(add.n)]]
    active[add.idx] <- TRUE

    # Reset acceleration state after active-set change.
    state <- .init_operator_state(prep$p, prep$k, intercept, B_init = B.cur)
  }

  if (ws_final_full) {
    core.full <- .run_operator_iterations(
      Y = prep$Y,
      groups = prep$groups,
      group.names = prep$group.names,
      XX = prep$XX,
      XY = prep$XY,
      X.list = prep$X.list,
      p = prep$p,
      k = prep$k,
      samp.sizes = prep$samp.sizes,
      lambda = lambda,
      gamma = gamma,
      mu = mu,
      tol = tol,
      num.it = ws_final_it,
      lam.max = lam.max,
      intercept = intercept,
      penalty.factors = prep$penalty.factors,
      edge.u = edge.u,
      edge.v = edge.v,
      edge.w = edge.w,
      edge.block = edge.block,
      state = .init_operator_state(prep$p, prep$k, intercept, B_init = state$B.old),
      trace_state = trace_state,
      stage = "ws_full",
      prep = prep
    )
    total.iter <- total.iter + core.full$iter
    state <- core.full$state
    converged <- converged || core.full$converged
    active[] <- TRUE
  }

  list(
    B = state$B.old,
    iterations = total.iter,
    converged = converged,
    active_edges = sum(active),
    total_edges = num.edges,
    screening_info = screen_info
  )
}

#' Fused L1 solver (operator + working-set edge generation).
#'
#' @param X matrix of covariates (n by p)
#' @param Y response vector (length n)
#' @param groups group indicators (length n)
#' @param lambda sparsity hyperparameter (scalar)
#' @param gamma fusion hyperparameter (scalar)
#' @param G pairwise fusion weight matrix (k by k)
#' @param mu smoothness parameter for proximal updates
#' @param tol optimization tolerance
#' @param num.it maximum iterations
#' @param lam.max optional precomputed maximal eigenvalue
#' @param c.flag Kept for API parity; ignored in this implementation
#' @param intercept whether to include per-group intercept terms
#' @param penalty.factors per-feature sparsity weights
#' @param conserve.memory whether to avoid materializing XX/XY lists
#' @param scaling whether to scale losses by group size
#' @param edge.block number of edges processed per block
#' @param ws_init_edges initial number of active edges
#' @param ws_add_edges max number of edges added per outer round
#' @param ws_max_outer max number of working-set outer rounds
#' @param ws_inner_it inner-iteration budget per outer round
#' @param ws_violation_tol edge-violation threshold for activation
#' @param ws_final_full if TRUE, run a final full-edge refinement pass
#' @param ws_final_it iteration budget for final full-edge refinement
#' @param screening edge screening strategy: "none" or "grad_zero"
#' @param screen_margin safety margin in [0,1] for `screening = "grad_zero"`
#' @param screen_max_drop_frac maximal fraction of edges allowed to be screened
#' @param screen_min_keep minimum number of edges to keep after screening
#' @param diagnostics whether to record per-iteration objective/KKT diagnostics
#' @param trace_every record diagnostics every `trace_every` local iterations
#'
#' @return Coefficient matrix (p by k, or (p+1) by k if intercept=TRUE)
#' @export
fusedLassoProximalNewWorkingSet <- function(
  X,
  Y,
  groups,
  lambda,
  gamma,
  G,
  mu = 1e-04,
  tol = 1e-06,
  num.it = 1000,
  lam.max = NULL,
  c.flag = FALSE,
  intercept = TRUE,
  penalty.factors = NULL,
  conserve.memory = NULL,
  scaling = TRUE,
  edge.block = 256L,
  ws_init_edges = 256L,
  ws_add_edges = 256L,
  ws_max_outer = 8L,
  ws_inner_it = NULL,
  ws_violation_tol = 1e-06,
  ws_final_full = TRUE,
  ws_final_it = NULL,
  screening = c("none", "grad_zero"),
  screen_margin = 0,
  screen_max_drop_frac = 1,
  screen_min_keep = 0L,
  diagnostics = FALSE,
  trace_every = 1L
) {
  .validate_l1new_common(lambda, gamma, mu, edge.block)
  edge.block <- as.integer(edge.block)
  screening <- match.arg(screening)
  if (!is.numeric(screen_margin) || length(screen_margin) != 1L || is.na(screen_margin) ||
      screen_margin < 0 || screen_margin > 1) {
    stop("Parameter 'screen_margin' must be a scalar in [0, 1].")
  }
  if (!is.numeric(screen_max_drop_frac) || length(screen_max_drop_frac) != 1L ||
      is.na(screen_max_drop_frac) || screen_max_drop_frac < 0 || screen_max_drop_frac > 1) {
    stop("Parameter 'screen_max_drop_frac' must be a scalar in [0, 1].")
  }
  if (!is.numeric(screen_min_keep) || length(screen_min_keep) != 1L ||
      is.na(screen_min_keep) || screen_min_keep < 0) {
    stop("Parameter 'screen_min_keep' must be a non-negative integer.")
  }
  screen_min_keep <- as.integer(screen_min_keep)

  if (c.flag) {
    warning("fusedLassoProximalNewWorkingSet does not use c.flag; running R implementation.")
  }
  if (is.null(conserve.memory)) {
    conserve.memory <- (ncol(X) >= 10000)
  }
  if (!is.logical(conserve.memory) || length(conserve.memory) != 1L) {
    stop("Parameter 'conserve.memory' must be TRUE/FALSE or NULL.")
  }

  prep <- .prepare_l1new_inputs(
    X = X,
    Y = Y,
    groups = groups,
    G = G,
    intercept = intercept,
    penalty.factors = penalty.factors,
    conserve.memory = conserve.memory,
    scaling = scaling
  )

  trace_state <- .init_trace_state(diagnostics = diagnostics, trace_every = trace_every)

  fit.raw <- .genFusedLassoProximalOperatorWS(
    prep = prep,
    lambda = lambda,
    gamma = gamma,
    mu = mu,
    tol = tol,
    num.it = num.it,
    lam.max = lam.max,
    intercept = intercept,
    edge.block = edge.block,
    ws_init_edges = ws_init_edges,
    ws_add_edges = ws_add_edges,
    ws_max_outer = ws_max_outer,
    ws_inner_it = ws_inner_it,
    ws_violation_tol = ws_violation_tol,
    ws_final_full = ws_final_full,
    ws_final_it = ws_final_it,
    screening = screening,
    screen_margin = screen_margin,
    screen_max_drop_frac = screen_max_drop_frac,
    screen_min_keep = screen_min_keep,
    trace_state = trace_state
  )

  .finalize_l1new_fit(
    fit.raw = fit.raw,
    prep = prep,
    lambda = lambda,
    intercept = intercept,
    trace_state = trace_state,
    screening_info = fit.raw$screening_info
  )
}
