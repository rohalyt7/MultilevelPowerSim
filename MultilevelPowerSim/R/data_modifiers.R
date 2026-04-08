

#' Set Outcome Distribution
#'
#' @description Replaces the existing outcome variable with values drawn from 
#' a specified probability distribution.
#'
#' @param data Dataframe. The simulated dataset.
#' @param outcome Character. The name of the column to replace.
#' @param dist Character. "normal", "binary", "poisson", "gamma", or "uniform".
#' @param ... Additional arguments passed to the distribution functions (e.g., lambda, shape).
#'
#' @return The dataframe with the modified outcome variable.
#' @export
set_outcome_distribution <- function(data, outcome = "Y", dist = "normal", ...) {
  n <- nrow(data)
  args <- list(...)
  
  # Helper to handle default values without needing rlang/%||%
  get_arg <- function(name, default) {
    if (name %in% names(args)) return(args[[name]]) else return(default)
  }
  
  data[[outcome]] <- switch(
    dist,
    "normal"  = rnorm(n, mean = get_arg("mean", 0), sd = get_arg("sd", 1)),
    "binary"  = rbinom(n, 1, prob = get_arg("prob", 0.5)),
    "poisson" = rpois(n, lambda = get_arg("lambda", 1)),
    "gamma"   = rgamma(n, shape = get_arg("shape", 2), scale = get_arg("scale", 1)),
    "uniform" = runif(n, min = get_arg("min", 0), max = get_arg("max", 1)),
    stop("Unsupported distribution type.")
  )
  
  message("Replaced outcome with ", dist, " distribution.")
  return(data)
}


#' Add Level-2 Covariate
#'
#' @description Generates a covariate that varies only at the cluster level and 
#' optionally injects its effect into the outcome Y.
#'
#' @param data Dataframe.
#' @param id_var Character. The column identifying clusters (e.g., "participant").
#' @param cov_name Character. Name for the new covariate.
#' @param dist Character. "normal", "uniform", or "binary".
#' @param mean Numeric. Mean for normal/binary, or min for uniform.
#' @param sd Numeric. SD for normal, or max for uniform.
#' @param beta Numeric. The effect size to add to the outcome Y.
#' @param outcome Character. The outcome column to modify.
#'
#' @export
add_level2_covariate <- function(data, id_var, cov_name = "Z",
                                 dist = "normal", mean = 0, sd = 1,
                                 beta = 0, outcome = NULL) {
  
  ids <- unique(data[[id_var]])
  n_ids <- length(ids)
  
  Z <- if (dist == "normal") {
    rnorm(n_ids, mean = mean, sd = sd)
  } else if (dist == "uniform") {
    runif(n_ids, min = mean, max = sd)
  } else if (dist == "binary") {
    rbinom(n_ids, 1, mean)
  } else {
    stop("Unsupported distribution.")
  }
  
  cov_df <- data.frame(id = ids, temp_Z = Z)
  names(cov_df) <- c(id_var, cov_name)
  
  data <- merge(data, cov_df, by = id_var)
  
  if (!is.null(outcome) && beta != 0) {
    data[[outcome]] <- data[[outcome]] + beta * data[[cov_name]]
  }
  
  return(data)
}


#' Add Bounded Variable
#'
#' @description Adds a variable (like time) that is bounded between two values.
#'
#' @param data Dataframe.
#' @param id_var Character. Cluster ID variable.
#' @param var_name Character. Name for the new variable.
#' @param lower Numeric. Lower bound.
#' @param upper Numeric. Upper bound.
#' @param n_points Integer. If provided, overrides cluster size.
#' @param random Logical. If TRUE, samples randomly; if FALSE, uses seq().
#' @param center Logical. If TRUE, centers the variable.
#' @param scale Logical. If TRUE, scales to unit variance.
#' @param outcome Character. Outcome to modify.
#' @param beta Numeric. Effect size.
#'
#' @export
add_bounded_variable <- function(data, id_var, var_name = "time",
                                 lower = 0, upper = 10, n_points = NULL,
                                 random = FALSE, center = FALSE, scale = FALSE,
                                 outcome = NULL, beta = 0) {
  
  # Split-Apply-Combine pattern
  data_list <- lapply(split(data, data[[id_var]]), function(cluster_data) {
    n_obs <- nrow(cluster_data)
    n <- if (!is.null(n_points)) n_points else n_obs
    
    t_vals <- if (!random) {
      seq(lower, upper, length.out = n)
    } else {
      sort(runif(n, min = lower, max = upper))
    }
    
    if (center) t_vals <- t_vals - mean(t_vals)
    if (scale) t_vals <- as.numeric(scale(t_vals, center = FALSE))
    
    cluster_data[[var_name]] <- t_vals
    
    if (!is.null(outcome) && beta != 0) {
      cluster_data[[outcome]] <- cluster_data[[outcome]] + beta * cluster_data[[var_name]]
    }
    return(cluster_data)
  })
  
  return(do.call(rbind, data_list))
}


#' Add Correlation to Predictors
#'
#' @description Uses a multivariate normal approach to impose a specific 
#' correlation structure on existing numeric predictors.
#'
#' @param df Dataframe.
#' @param predictors Vector. Names of columns to correlate.
#' @param cor_matrix Matrix. A valid correlation matrix.
#' @param betas Named vector. If provided, regenerates Y to maintain effect sizes.
#' @param group_var Character. Categorical variable to handle in Y regeneration.
#'
#' @export
add_correlation <- function(df, predictors, cor_matrix, betas = NULL, group_var = "Group") {
  
  if (length(predictors) < 2) {
    warning("add_correlation() needs at least two predictors.")
    return(df)
  }
  
  # Extract numeric data
  X <- as.matrix(df[, predictors, drop = FALSE])
  means <- colMeans(X)
  sds <- apply(X, 2, sd)
  
  # Reconstruct covariance matrix
  cov_matrix <- diag(sds) %*% cor_matrix %*% diag(sds)
  
  # Generate correlated data
  X_corr <- MASS::mvrnorm(n = nrow(df), mu = means, Sigma = cov_matrix)
  df[, predictors] <- X_corr
  
  # Optional: Regenerate Y because changing X changes the relationship with Y
  if (!is.null(betas)) {
    intercept <- if ("(Intercept)" %in% names(betas)) betas["(Intercept)"] else 0
    eta <- rep(intercept, nrow(df))
    
    for (pred in predictors) {
      if (pred %in% names(betas)) eta <- eta + betas[pred] * df[[pred]]
    }
    
    # Handle categorical groups
    if (group_var %in% names(df)) {
      lvls <- unique(df[[group_var]])
      for (lvl in lvls[-1]) {
        term <- paste0(group_var, lvl)
        if (term %in% names(betas)) {
          eta <- eta + betas[term] * (df[[group_var]] == lvl)
        }
      }
    }
    df$Y <- eta + rnorm(nrow(df), 0, 1)
  }
  
  return(df)
}


#' Add Heteroskedasticity
#'
#' @description Modifies Y such that the residual variance is a function 
#' of a specific predictor (variance increases/decreases with X).
#'
#' @param data Dataframe.
#' @param eta Numeric vector. The linear predictor (X %*% beta).
#' @param predictor Character. The name of the predictor driving the variance.
#' @param sigma0 Numeric. Baseline standard deviation.
#' @param factor Numeric. Magnitude of the variance increase.
#'
#' @export
add_heteroskedasticity <- function(data, eta, predictor, sigma0 = 1, factor = 0.5) {
  # sigma_i = baseline * (1 + factor * X)
  sigma_i <- sigma0 * (1 + factor * data[[predictor]])
  data$Y <- eta + rnorm(nrow(data), mean = 0, sd = sigma_i)
  
  return(data)
}


#' Scale All Numeric Predictors
#' @export
scale_predictors <- function(data, outcome = "Y") {
  numeric_cols <- setdiff(names(data)[sapply(data, is.numeric)], outcome)
  data[numeric_cols] <- lapply(data[numeric_cols], scale)
  return(data)
}

#' Advanced Centering and Scaling
#' @description Supports group-mean and grand-mean centering.
#' @export
center_scale_vars <- function(data, id_var, group_mean = NULL, 
                              grand_mean = NULL, scale = TRUE, suffix = "_cs") {
  
  df <- as.data.frame(data)
  
  # Group-mean centering (Within-cluster)
  if (!is.null(group_mean)) {
    for (var in group_mean) {
      gm_name <- paste0(var, suffix)
      x <- as.numeric(df[[var]])
      g <- df[[id_var]]
      df[[gm_name]] <- ave(x, g, FUN = function(xi) xi - mean(xi, na.rm = TRUE))
      if (scale) df[[gm_name]] <- as.numeric(scale(df[[gm_name]], center = FALSE))
    }
  }
  
  # Grand-mean centering
  if (!is.null(grand_mean)) {
    for (var in grand_mean) {
      gm_name <- paste0(var, suffix)
      df[[gm_name]] <- as.numeric(df[[var]] - mean(df[[var]], na.rm = TRUE))
      if (scale) df[[gm_name]] <- as.numeric(scale(df[[gm_name]], center = FALSE))
    }
  }
  
  return(df)
}