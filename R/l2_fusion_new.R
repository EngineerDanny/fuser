# Experimental L2 fusion implementation focused on memory efficiency.
# Changes vs current implementation:
# 1) Build fusion constraints only for non-zero edges in G.
# 2) Construct sparse matrices from triplets (i, j, x) in one shot.

#' Generate transformed matrices for fused L2 using sparse edge construction.
#'
#' @param X Covariates matrix (n by p).
#' @param Y Response vector (length n).
#' @param groups Group indicators (length n).
#' @param G Fusion strength matrix (K by K).
#' @param intercept Whether to include a per-group intercept.
#' @param penalty.factors Penalty factors for covariates.
#' @param scaling Whether to scale each subgroup by sqrt(n_k).
#' @param include.fusion Whether to construct fusion constraints from G.
#'
#' @return A list with X, Y, X.fused, penalty, edges, and group.names.
generateBlockDiagonalMatricesNew <- function(
  X,
  Y,
  groups,
  G,
  intercept = FALSE,
  penalty.factors = rep(1, dim(X)[2]),
  scaling = FALSE,
  include.fusion = TRUE
) {
  group.names <- sort(unique(groups))
  num.groups <- length(group.names)

  if (!is.matrix(G) || nrow(G) != num.groups || ncol(G) != num.groups) {
    stop("G must be a square matrix with dimensions equal to the number of groups.")
  }

  if (intercept) {
    X <- cbind(X, matrix(1, nrow(X), 1))
  }

  n <- nrow(X)
  p.eff <- ncol(X)

  # Keep only active fusion edges (upper triangle, non-zero tau) when requested.
  if (include.fusion) {
    edges <- which(upper.tri(G) & (G != 0), arr.ind = TRUE)
  } else {
    edges <- matrix(integer(0), nrow = 0, ncol = 2)
  }
  num.edges <- nrow(edges)

  # Response vector is original data rows + optional fusion rows (zeros).
  new.y <- numeric(length(Y) + num.edges * p.eff)

  # Build block-diagonal X via triplets.
  i.list <- vector("list", num.groups)
  j.list <- vector("list", num.groups)
  x.list <- vector("list", num.groups)

  row.start <- 1L

  for (g.idx in seq_len(num.groups)) {
    group.inds <- which(groups == group.names[g.idx])
    n.g <- length(group.inds)

    row.range <- row.start:(row.start + n.g - 1L)
    col.base <- (g.idx - 1L) * p.eff

    block <- X[group.inds, , drop = FALSE]
    if (scaling) {
      scale.factor <- sqrt(n.g)
      block <- block / scale.factor
      new.y[row.range] <- Y[group.inds] / scale.factor
    } else {
      new.y[row.range] <- Y[group.inds]
    }

    # Dense block converted to triplets; zeros are filtered out.
    block.vec <- c(t(block))
    ii <- rep(row.range, each = p.eff)
    jj <- rep(col.base + seq_len(p.eff), times = n.g)
    keep <- (block.vec != 0)

    i.list[[g.idx]] <- ii[keep]
    j.list[[g.idx]] <- jj[keep]
    x.list[[g.idx]] <- block.vec[keep]

    row.start <- row.start + n.g
  }

  new.x <- Matrix::sparseMatrix(
    i = unlist(i.list, use.names = FALSE),
    j = unlist(j.list, use.names = FALSE),
    x = unlist(x.list, use.names = FALSE),
    dims = c(n, p.eff * num.groups)
  )

  # Build fusion matrix X.fused via triplets (only active edges).
  new.x.f <- Matrix::Matrix(0, num.edges * p.eff, p.eff * num.groups, sparse = TRUE)

  if (num.edges > 0) {
    feat.idx <- seq_len(p.eff)
    if (intercept) {
      feat.idx <- feat.idx[-length(feat.idx)] # Do not fuse intercept.
    }

    if (length(feat.idx) > 0) {
      i.f.list <- vector("list", num.edges)
      j.f.list <- vector("list", num.edges)
      x.f.list <- vector("list", num.edges)

      for (e.idx in seq_len(num.edges)) {
        g.i <- edges[e.idx, 1]
        g.j <- edges[e.idx, 2]
        tau <- G[g.i, g.j]
        coeff <- sqrt(tau)

        rows <- (e.idx - 1L) * p.eff + feat.idx
        cols.i <- (g.i - 1L) * p.eff + feat.idx
        cols.j <- (g.j - 1L) * p.eff + feat.idx

        i.f.list[[e.idx]] <- c(rows, rows)
        j.f.list[[e.idx]] <- c(cols.i, cols.j)
        x.f.list[[e.idx]] <- c(rep(coeff, length(feat.idx)), rep(-coeff, length(feat.idx)))
      }

      new.x.f <- Matrix::sparseMatrix(
        i = unlist(i.f.list, use.names = FALSE),
        j = unlist(j.f.list, use.names = FALSE),
        x = unlist(x.f.list, use.names = FALSE),
        dims = c(num.edges * p.eff, p.eff * num.groups)
      )
    }
  }

  if (intercept) {
    penalty <- rep(c(penalty.factors, 0), num.groups)
  } else {
    penalty <- rep(penalty.factors, num.groups)
  }

  return(list(
    X = new.x,
    Y = Matrix::Matrix(new.y, ncol = 1, sparse = TRUE),
    X.fused = new.x.f,
    penalty = penalty,
    edges = edges,
    group.names = group.names
  ))
}


#' Optimise fused L2 with active-edge transformed matrices.
#'
#' @param X Covariates matrix.
#' @param y Response vector.
#' @param groups Group indicators.
#' @param lambda Sparsity hyperparameter.
#' @param G Fusion matrix.
#' @param gamma Fusion multiplier.
#' @param scaling Whether to scale each subgroup.
#' @param ... Additional parameters passed to glmnet.
#'
#' @return Coefficient matrix (p by k).
fusedL2DescentGLMNetNew <- function(
  X,
  y,
  groups,
  lambda = NULL,
  G = NULL,
  gamma = 1,
  scaling = FALSE,
  ...
) {
  if (!is.numeric(gamma) || length(gamma) != 1 || gamma < 0) {
    stop("gamma must be a non-negative scalar.")
  }

  if (is.null(G)) {
    k <- length(sort(unique(groups)))
    G <- matrix(1, k, k)
  }

  include.fusion <- (gamma > 0)
  transformed.data <- generateBlockDiagonalMatricesNew(
    X, y, groups, G,
    scaling = scaling,
    include.fusion = include.fusion
  )
  transformed.x <- transformed.data$X
  transformed.x.f <- transformed.data$X.fused
  transformed.y <- transformed.data$Y

  if (include.fusion && nrow(transformed.x.f) > 0) {
    transformed.x.f <- transformed.x.f * sqrt(gamma * (nrow(transformed.x) + nrow(transformed.x.f)))
    transformed.x <- rbind(transformed.x, transformed.x.f)
  }

  group.names <- sort(unique(groups))
  num.groups <- length(group.names)

  if (scaling) {
    correction.factor <- num.groups / nrow(transformed.x)
  } else {
    correction.factor <- length(groups) / nrow(transformed.x)
  }

  glmnet.result <- glmnet::glmnet(
    transformed.x,
    transformed.y,
    standardize = FALSE,
    intercept = FALSE,
    lambda = lambda * correction.factor,
    ...
  )

  coef.temp <- coef(glmnet.result, s = lambda * correction.factor)
  beta.mat <- matrix(
    coef.temp[2:length(coef.temp)],
    nrow = ncol(transformed.x) / num.groups,
    ncol = num.groups
  )
  colnames(beta.mat) <- group.names

  return(beta.mat)
}
