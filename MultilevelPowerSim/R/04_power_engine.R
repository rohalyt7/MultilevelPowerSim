# Dependency Check --------------------------------------------------------
if (!requireNamespace("lme4", quietly = TRUE)) stop("Package 'lme4' is required.")
if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

# --- Internal Helpers ---

#' Extract Model Diagnostics
#'
#' @description Pulls convergence, singularity, and point estimates from an lme4 model.
#'
#' @param fit An lme4 model object (lmer or glmer).
#' @param effect_name Character. The name of the fixed effect to extract.
#'
#' @return A list containing singularity status, convergence messages, and the estimate.
#' @keywords internal
extract_diagnostics <- function(fit, effect_name) {
  list(
    is_singular = lme4::isSingular(fit),
    has_converge_warn = !is.null(fit@optinfo$conv$lme4$messages),
    est = lme4::fixef(fit)[effect_name]
  )
}

#' Run Individual Simulation Iterations
#'
#' @description The workhorse loop that handles data generation, model fitting, 
#' and captures warnings/errors for a single sample size (n).
#'
#' @param n Integer. Current sample size.
#' @param n_sim Integer. Total number of simulations to run.
#' @param data_generator Function. A user-defined function that accepts 'n' and returns data.
#' @param fit_function Function. A user-defined function that fits a model to the data.
#' @param effect_name Character. The coefficient name to be tested for power.
#' @param alpha Numeric. The significance threshold (usually 0.05).
#'
#' @return A list containing summary statistics (power, fail rate, etc.) and raw diagnostics.
#' @keywords internal
run_simulations <- function(n, n_sim, data_generator, fit_function, effect_name, alpha) {
  pvals <- numeric(n_sim)
  ests <- numeric(n_sim)
  singular_flags <- logical(n_sim)
  converge_flags <- logical(n_sim)
  fail_count <- 0
  messages <- character()
  
  for (sim in seq_len(n_sim)) {
    # 1. Data Generation (Wrapped to catch logic errors in the generator)
    dat <- tryCatch(data_generator(n), error = function(e) { 
      messages <<- c(messages, paste("Data Error:", e$message))
      return(NULL) 
    })
    if (is.null(dat)) { fail_count <- fail_count + 1; next }
    
    # 2. Model Fitting (Captures warnings for singularity and convergence)
    fit <- tryCatch({
      withCallingHandlers(
        fit_function(dat),
        warning = function(w) {
          msg <- conditionMessage(w)
          if (grepl("converg", msg, ignore.case = TRUE)) converge_flags[sim] <<- TRUE
          if (grepl("singular|boundary", msg, ignore.case = TRUE)) singular_flags[sim] <<- TRUE
          invokeRestart("muffleWarning")
        }
      )
    }, error = function(e) { 
      messages <<- c(messages, paste("Model Error:", e$message))
      return(NULL) 
    })
    
    if (is.null(fit)) { fail_count <- fail_count + 1; next }
    
    # 3. Parameter Extraction
    fe <- lme4::fixef(fit)
    ests[sim] <- if (effect_name %in% names(fe)) fe[effect_name] else NA_real_
    
    # 4. P-Value Extraction (Handles varying column names like Pr(>|t|) or Pr(>|z|))
    pval <- tryCatch({
      coefs <- summary(fit)$coefficients
      row_idx <- grep(paste0("^", effect_name, "$"), rownames(coefs))
      if (length(row_idx) == 1) {
        p_col <- which(colnames(coefs) %in% c("Pr(>|t|)", "Pr(>|z|)"))
        if (length(p_col) > 0) coefs[row_idx, p_col] else NA_real_
      } else { NA_real_ }
    }, error = function(e) NA_real_)
    
    pvals[sim] <- pval
  }
  
  summary_stats <- c(
    power = mean(pvals < alpha, na.rm = TRUE),
    failed = fail_count / n_sim,
    pct_singular = mean(singular_flags, na.rm = TRUE),
    pct_converge_warn = mean(converge_flags, na.rm = TRUE),
    mean_est = mean(ests, na.rm = TRUE),
    sd_est = sd(ests, na.rm = TRUE)
  )
  
  diagnostics <- data.frame(
    n = n, sim = seq_len(n_sim), est = ests, pval = pvals,
    singular = singular_flags, converge_warn = converge_flags
  )
  
  return(list(summary = summary_stats, diagnostics = diagnostics, messages = messages))
}

# --- Main Functions ---

#' Main Power Analysis Driver
#'
#' @description Orchestrates the simulation-based power analysis across a 
#' vector of sample sizes.
#'
#' @param data_generator Function. Accepts 'n' and returns a simulated dataframe.
#' @param n_values Numeric vector. The sample sizes to evaluate.
#' @param fit_function Function. Accepts data and returns a fitted model.
#' @param effect_name Character. The name of the coefficient to calculate power for.
#' @param true_beta Numeric. Optional 'true' effect size for calculating bias/attenuation.
#' @param n_sim Integer. Number of simulations to perform at each sample size.
#' @param alpha Numeric. Alpha level for significance (default 0.05).
#' @param plot Logical. Should a power curve be automatically generated?
#' @param verbose Logical. Should progress messages be printed to the console?
#'
#' @return A list containing the 'power_table', 'power_plot', 'diagnostics', and 'failure_log'.
#' @export
power_analysis <- function(data_generator, n_values, fit_function,
                           effect_name, true_beta = NULL,
                           n_sim = 200, alpha = 0.05,
                           plot = TRUE, verbose = TRUE) {
  
  results <- data.frame(
    n = n_values, power = NA_real_, failed = NA_real_,
    pct_singular = NA_real_, pct_converge_warn = NA_real_,
    mean_est = NA_real_, sd_est = NA_real_, attenuation = NA_real_
  )
  
  failure_log <- vector("list", length(n_values))
  diagnostics_list <- vector("list", length(n_values))
  
  for (i in seq_along(n_values)) {
    n <- n_values[i]
    if (verbose) message("Running power analysis for n = ", n, "...")
    
    sim_out <- run_simulations(n, n_sim, data_generator, fit_function, effect_name, alpha)
    
    stat_names <- names(sim_out$summary)
    results[i, stat_names] <- sim_out$summary
    
    if (!is.null(true_beta)) results$attenuation[i] <- sim_out$summary["mean_est"] / true_beta
    
    failure_log[[i]] <- sim_out$messages
    diagnostics_list[[i]] <- sim_out$diagnostics
  }
  
  p_plot <- if (plot) plot_power_curve(results) else NULL
  
  return(list(
    power_table = results,
    power_plot = p_plot,
    diagnostics = diagnostics_list,
    failure_log = failure_log
  ))
}

#' Plot Power Curve
#'
#' @description Visualizes the relationship between sample size and statistical power.
#'
#' @param results Dataframe. The 'power_table' from 'power_analysis()'.
#'
#' @return A ggplot2 object.
#' @export
plot_power_curve <- function(results) {
  ggplot2::ggplot(results, ggplot2::aes(x = n, y = power)) +
    ggplot2::geom_line(linewidth = 1, color = "#2E86AB") +
    ggplot2::geom_point(size = 2, color = "#2E86AB") +
    ggplot2::labs(title = "Simulation-Based Power Curve",
                  x = "Sample Size (n)", y = "Estimated Power") +
    ggplot2::scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.1)) +
    ggplot2::theme_classic(base_size = 13)
}

#' Combine Power Tables
#'
#' @description Merges multiple power analysis results for comparison 
#' (e.g., across different ICC levels or effect sizes).
#'
#' @param power_tables List. A list of dataframes from 'power_analysis()'.
#' @param effect_sizes Numeric vector. The effect sizes corresponding to each table.
#' @param icc_values Numeric vector. The ICC values corresponding to each table.
#'
#' @return A single combined dataframe.
#' @export
combine_power_tables <- function(power_tables, effect_sizes, icc_values) {
  combined_list <- Map(function(tbl, eff, icc) {
    tbl$Effect.Size <- eff
    tbl$ICC <- icc
    tbl
  }, power_tables, effect_sizes, icc_values)
  
  do.call(rbind, combined_list)
}

#' Plot Power Heatmap
#'
#' @description Creates a tile plot of power across two dimensions 
#' (e.g., Effect Size and ICC).
#'
#' @param df Dataframe. A combined power table (see 'combine_power_tables').
#' @param palette Character. A viridis color palette name (default "plasma").
#'
#' @return A ggplot2 object.
#' @export
plot_power_heatmap <- function(df, palette = "plasma") {
  ggplot2::ggplot(df, ggplot2::aes(x = factor(Effect.Size), y = factor(ICC), fill = power)) +
    ggplot2::geom_tile(color = "white", linewidth = 0.5) +
    ggplot2::scale_fill_viridis_c(option = palette, direction = 1, name = "Power") +
    ggplot2::labs(title = "Power by Effect Size and ICC",
                  x = "Effect Size", y = "ICC Level") +
    ggplot2::theme_minimal(base_size = 13)
}
