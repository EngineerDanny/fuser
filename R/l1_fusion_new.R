# Backward-compatible wrapper over the split L1-new solvers.
# - Operator-only path: fusedLassoProximalNewOperator()
# - Working-set path:   fusedLassoProximalNewWorkingSet()

#' Fused L1 solver (compatibility wrapper).
#'
#' @param ... Arguments forwarded to the selected solver.
#' @param working_set If TRUE, use working-set solver, otherwise operator-only solver.
#' @param ws_init_edges Initial number of active edges for working-set mode.
#' @param ws_add_edges Max number of edges added per outer working-set round.
#' @param ws_max_outer Max number of outer working-set rounds.
#' @param ws_inner_it Inner-iteration budget per outer round.
#' @param ws_violation_tol Violation threshold for edge activation.
#' @param ws_final_full If TRUE, run a final full-edge refinement pass.
#' @param ws_final_it Iteration budget for final full-edge refinement.
#' @param screening Edge screening strategy in working-set mode: "none" or "grad_zero".
#' @param screen_margin Screening margin in [0,1] for `screening = "grad_zero"`.
#' @param screen_max_drop_frac Maximum fraction of edges allowed to be screened.
#' @param screen_min_keep Minimum number of edges to keep after screening.
#'
#' @return Coefficient matrix.
#' @export
fusedLassoProximalNew <- function(
  ...,
  working_set = FALSE,
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
  screen_min_keep = 0L
) {
  if (!is.logical(working_set) || length(working_set) != 1L) {
    stop("Parameter 'working_set' must be TRUE/FALSE.")
  }
  screening <- match.arg(screening)

  if (working_set) {
    fusedLassoProximalNewWorkingSet(
      ...,
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
      screen_min_keep = screen_min_keep
    )
  } else {
    fusedLassoProximalNewOperator(...)
  }
}
