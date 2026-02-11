.genFusedLassoProximalOperatorOnly <- function(
  prep,
  lambda,
  gamma,
  mu,
  tol,
  num.it,
  lam.max,
  intercept,
  edge.block,
  trace_state = NULL
) {
  if (is.null(lam.max)) {
    lam.max <- .compute_lam_max(prep$XX, prep$X.list)
  }

  init.state <- .init_operator_state(prep$p, prep$k, intercept)
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
    num.it = num.it,
    lam.max = lam.max,
    intercept = intercept,
    penalty.factors = prep$penalty.factors,
      edge.u = prep$edge.u,
      edge.v = prep$edge.v,
      edge.w = prep$edge.w,
      edge.block = edge.block,
      state = init.state,
      trace_state = trace_state,
      stage = "operator",
      prep = prep
    )

  list(
    B = core$B,
    iterations = core$iter,
    converged = core$converged,
    active_edges = length(prep$edge.u),
    total_edges = length(prep$edge.u)
  )
}

#' Fused L1 solver (operator-only, no dense fusion matrix materialization).
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
#' @param diagnostics whether to record per-iteration objective/KKT diagnostics
#' @param trace_every record diagnostics every `trace_every` local iterations
#'
#' @return Coefficient matrix (p by k, or (p+1) by k if intercept=TRUE)
#' @export
fusedLassoProximalNewOperator <- function(
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
  diagnostics = FALSE,
  trace_every = 1L
) {
  .validate_l1new_common(lambda, gamma, mu, edge.block)
  edge.block <- as.integer(edge.block)

  if (c.flag) {
    warning("fusedLassoProximalNewOperator does not use c.flag; running R implementation.")
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

  fit.raw <- .genFusedLassoProximalOperatorOnly(
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
  )

  .finalize_l1new_fit(
    fit.raw = fit.raw,
    prep = prep,
    lambda = lambda,
    intercept = intercept,
    trace_state = trace_state
  )
}
