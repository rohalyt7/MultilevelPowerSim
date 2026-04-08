# Dependency Check --------------------------------------------------------

if (!requireNamespace("MASS", quietly = TRUE)) {
  stop("Package 'MASS' is required for multivariate normal simulation. Please install it.")
}

# Note: 'stats' is a base package, so requireNamespace isn't strictly needed,
# but we use functions like rnorm, rbinom, and model.matrix from it.

#' Simulate Basic Linear Data
#'
#' @description Generates a foundation dataset with various predictor types 
#' (Normal, Binary, Categorical) and a continuous outcome Y based on a 
#' linear combination of those predictors.
#'
#' @param n Integer. Number of observations.
#' @param betas Numeric vector. The "true" coefficients, including the Intercept.
#' @param predictors Named list. Each element defines a predictor (type, and parameters).
#' @param error_sd Numeric. Standard deviation of the Gaussian noise (residual error).
#'
#' @return A list containing the generated 'data', the 'model_matrix' (X), 
#' underlying 'eta' (linear predictor), and the 'betas' used.
#' @export
simulate_data <- function(n, betas, predictors, error_sd = 1) {
  # Step 1: Initialize dataframe
  df <- data.frame(row_id = seq_len(n))
  
  # Step 2: Generate predictors
  for (pred_name in names(predictors)) {
    spec <- predictors[[pred_name]]
    
    if (spec$type == "normal") {
      df[[pred_name]] <- rnorm(n, mean = spec$mean, sd = spec$sd)
    } else if (spec$type == "binary") {
      df[[pred_name]] <- rbinom(n, 1, spec$prob)
    } else if (spec$type == "categorical") {
      # Determine levels from names or indices
      levels_vec <- if (!is.null(names(spec$probs))) names(spec$probs) else seq_along(spec$probs)
      df[[pred_name]] <- factor(
        sample(levels_vec, n, replace = TRUE, prob = spec$probs),
        levels = levels_vec
      )
    }
  }
  
  # Step 3: Construct model matrix (X)
  # Exclude row_id from the design matrix
  X <- model.matrix(~ . - row_id, data = df)
  
  # Validation: Check if beta length matches the design matrix columns
  if (length(betas) != ncol(X)) {
    stop("Length of betas (", length(betas), ") does not match ncol(X) (", ncol(X), ").\n",
         "Columns: ", paste(colnames(X), collapse = ", "))
  }
  
  # Step 4: Generate linear predictor (eta) and outcome (Y)
  eta <- as.numeric(X %*% betas) 
  y <- eta + rnorm(n, 0, error_sd) 
  df$Y <- y
  
  return(list(
    data = df,
    model_matrix = X,
    betas = betas,
    eta = eta,
    error_sd = error_sd  # Crucial for the ICC derivation in the next step
  ))
}


#' Simulate Clustered Data from Scratch
#'
#' @description Generates a multilevel dataset with random intercepts and 
#' random slopes for a specific predictor across clusters.
#'
#' @param n_clusters Integer. Number of clusters/participants.
#' @param obs_per_cluster Integer. Number of observations per cluster.
#' @param predictors Named list. Specification for predictors.
#' @param betas Numeric vector. Fixed effects coefficients.
#' @param slope_var Character. The name of the predictor that varies randomly.
#' @param error_sd Numeric. Level-1 residual standard deviation.
#' @param icc Numeric. Intraclass Correlation Coefficient (used to derive intercept variance).
#' @param tau1 Numeric. Standard deviation of the random slope.
#' @param cov01 Numeric. Covariance between random intercept and random slope.
#'
#' @export
simulate_clustered_data <- function(n_clusters,
                                    obs_per_cluster = 1,
                                    predictors,
                                    betas,
                                    slope_var,
                                    error_sd = 1,
                                    icc = 0.2,
                                    tau1 = 0.5,
                                    cov01 = 0) {
  
  if (!slope_var %in% names(predictors)) {
    stop("slope_var must be one of the predictors.")
  }
  
  # Step 1: Random effect covariance matrix (D)
  # Derive tau0 (intercept SD) from ICC
  tau0 <- sqrt((icc * error_sd^2) / (1 - icc))
  D <- matrix(c(tau0^2, cov01, cov01, tau1^2), nrow = 2)
  
  # Step 2: Generate random effects per cluster
  u <- MASS::mvrnorm(n_clusters, mu = c(0, 0), Sigma = D)
  colnames(u) <- c("u0", "u1")
  
  data_list <- vector("list", n_clusters)
  
  # Step 3: Loop through clusters
  for (i in 1:n_clusters) {
    cluster_df <- data.frame(
      participant = i,
      lapply(predictors, function(spec) {
        if (spec$type == "normal") {
          rnorm(obs_per_cluster, mean = spec$mean, sd = spec$sd)
        } else if (spec$type == "binary") {
          rbinom(obs_per_cluster, 1, spec$prob)
        } else if (spec$type == "categorical") {
          # Validation for categorical probabilities
          if (is.null(spec$probs)) stop("Categorical predictor must include 'probs'.")
          if (any(spec$probs < 0)) stop("Probabilities must be non-negative.")
          if (abs(sum(spec$probs) - 1) > 1e-6) stop("Probabilities must sum to 1.")
          
          levels_vec <- if (!is.null(names(spec$probs))) {
            names(spec$probs)
          } else {
            paste0("L", seq_along(spec$probs))
          }
          sample(levels_vec, obs_per_cluster, replace = TRUE, prob = spec$probs)
        }
      })
    )
    
    # Construct fixed effects model matrix
    X <- model.matrix(~ . - participant, data = cluster_df)
    
    if (length(betas) != ncol(X)) {
      stop("Length of betas does not match number of predictors (including intercept).")
    }
    
    # Calculate Y with fixed effects + random intercept + random slope
    eta <- as.numeric(X %*% betas)
    
    # Crucial: cast slope_var to numeric to avoid factor math errors
    slope_val <- as.numeric(cluster_df[[slope_var]])
    
    cluster_df$Y <- eta + u[i, "u0"] + u[i, "u1"] * slope_val +
      rnorm(obs_per_cluster, 0, error_sd)
    
    data_list[[i]] <- cluster_df
  }
  
  full_df <- do.call(rbind, data_list)
  
  return(list(
    data = full_df,
    betas = betas,
    u0 = u[, "u0"],
    u1 = u[, "u1"],
    D = D,
    error_sd = error_sd,
    icc = icc,
    n_clusters = n_clusters,
    obs_per_cluster = obs_per_cluster
  ))
}


#' Add Clustering and Random Slopes to an Existing Simulation
#'
#' @description Takes a simulation object (from simulate_data) and expands it 
#' by adding clusters (participants), random intercepts, and random slopes.
#'
#' @param sim List. The output object from 'simulate_data()'.
#' @param obs_per_participant Integer. Observations per cluster.
#' @param icc Numeric. Intraclass Correlation Coefficient.
#' @param tau1 Numeric. Standard deviation of the random slope.
#' @param cov01 Numeric. Covariance between random intercept and slope.
#' @param slope_var Character. Predictor name that carries the random slope.
#'
#' @export
add_cluster_random_slope <- function(sim,
                                     obs_per_participant = 1,
                                     icc,
                                     tau1,
                                     cov01 = 0,
                                     slope_var = NULL) {
  df <- sim$data
  n_clusters <- nrow(df)
  
  if (!slope_var %in% names(df)) {
    stop(paste0("slope_var '", slope_var, "' not found in data."))
  }
  
  # Derive tau0 from ICC
  tau0 <- sqrt((icc * sim$error_sd^2) / (1 - icc))
  D <- matrix(c(tau0^2, cov01, cov01, tau1^2), nrow = 2)
  
  # Generate random effects
  u <- MASS::mvrnorm(n_clusters, mu = c(0, 0), Sigma = D)
  colnames(u) <- c("u0", "u1")
  
  # Expand data rows for clustering
  df <- df[rep(1:n_clusters, each = obs_per_participant), , drop = FALSE]
  df$participant <- rep(1:n_clusters, each = obs_per_participant)
  
  # Update outcome Y: Original Y + u0 + (u1 * slope_var)
  # Ensure slope_var is numeric
  slope_val <- as.numeric(df[[slope_var]])
  df$Y <- df$Y + u[df$participant, "u0"] + u[df$participant, "u1"] * slope_val
  
  # Update the simulation object
  sim$data <- df
  sim$u0 <- u[, "u0"]
  sim$u1 <- u[, "u1"]
  sim$D <- D
  sim$icc <- icc
  sim$n_clusters <- n_clusters
  sim$obs_per_participant <- obs_per_participant
  
  return(sim)
}
