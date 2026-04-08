# Dependency Check --------------------------------------------------------
# This file primarily uses base stats functions.

#' Introduce Missing Data (MCAR, MAR, MNAR)
#'
#' @description Injects missing values (NAs) into numeric or categorical variables 
#' using standard missingness mechanisms. This is essential for testing model 
#' robustness against incomplete data.
#'
#' @param data Dataframe. The simulated dataset to modify.
#' @param vars Character vector. Names of numeric variables to make missing.
#' @param cat_vars Character vector. Names of categorical variables to make missing.
#' @param prop_missing Numeric vector. The proportion (0 to 1) of missingness for 
#' each variable. Must match the combined length of 'vars' and 'cat_vars'.
#' @param mechanism Character. The mechanism of missingness: "MCAR" (Random), 
#' "MAR" (Conditional on another variable), or "MNAR" (Conditional on itself).
#' @param ref_var Character. The name of the reference variable required for 
#' the MAR mechanism.
#'
#' @return A dataframe with NAs injected into the specified columns.
#' @export
add_missing_data <- function(data,
                             vars = NULL,
                             cat_vars = NULL,
                             prop_missing,
                             mechanism = c("MCAR", "MAR", "MNAR"),
                             ref_var = NULL) {
  
  mechanism <- match.arg(mechanism)
  total_vars <- c(vars, cat_vars)
  
  # Validation: Ensure mapping between variables and proportions is 1:1
  if (length(prop_missing) != length(total_vars)) {
    stop("Length of prop_missing must match total number of vars + cat_vars.")
  }
  
  # Validation: MAR requires a 'cause' variable to be specified
  if (mechanism == "MAR" && is.null(ref_var)) {
    stop("MAR mechanism requires a reference variable (ref_var) to drive missingness.")
  }
  
  df <- data
  
  # 1. Process Numeric Variables
  if (!is.null(vars)) {
    for (i in seq_along(vars)) {
      var <- vars[i]
      prop <- prop_missing[i]
      
      # Determine which rows become NA based on the selected mechanism
      mask <- switch(
        mechanism,
        "MCAR" = mcar_mask(nrow(df), prop),
        "MAR"  = mar_mask_numeric(df[[ref_var]], prop),
        "MNAR" = mar_mask_numeric(df[[var]], prop)
      )
      df[[var]][mask] <- NA
    }
  }
  
  # 2. Process Categorical Variables
  if (!is.null(cat_vars)) {
    for (j in seq_along(cat_vars)) {
      var <- cat_vars[j]
      prop <- prop_missing[length(vars) + j]
      
      mask <- switch(
        mechanism,
        "MCAR" = mcar_mask(nrow(df), prop),
        "MAR"  = mar_mask_categorical(df[[ref_var]], prop),
        "MNAR" = mar_mask_categorical(df[[var]], prop)
      )
      
      # Handle factors safely: convert to character to inject NA, then rebuild factor
      # This prevents issues with "level" mismatches during NA injection
      tmp <- as.character(df[[var]])
      orig_levels <- levels(df[[var]])
      tmp[mask] <- NA
      df[[var]] <- factor(tmp, levels = orig_levels)
    }
  }
  
  return(df)
}

# --- Helper Functions (Internal) ---

#' MCAR Mask
#' @description Generates a random boolean mask where TRUE indicates missing.
#' @param n Integer. Number of observations.
#' @param prop Numeric. Target proportion of missingness.
#' @keywords internal
mcar_mask <- function(n, prop) {
  sample(c(TRUE, FALSE), size = n, replace = TRUE, prob = c(prop, 1 - prop))
}

#' MAR/MNAR Mask for Numeric Predictors
#' @description Uses a logistic/quantile approach to select the 'highest' values for missingness.
#' @param x Numeric vector. The variable driving the missingness.
#' @param prop Numeric. Target proportion of missingness.
#' @keywords internal
mar_mask_numeric <- function(x, prop) {
  # Handle edge case: if variable has no variance, default to random
  if (sd(x, na.rm = TRUE) == 0) return(mcar_mask(length(x), prop))
  
  # Transform values to a probability scale
  p <- stats::plogis(scale(x)) 
  # Create a deterministic threshold based on the desired proportion
  thresh <- stats::quantile(p, probs = 1 - prop, na.rm = TRUE)
  mask <- p > thresh
  return(as.logical(mask))
}

#' MAR/MNAR Mask for Categorical Predictors
#' @description Assigns missingness likelihood based on factor levels.
#' @param x Factor/Vector. The categorical variable.
#' @param prop Numeric. Target proportion of missingness.
#' @keywords internal
mar_mask_categorical <- function(x, prop) {
  x <- factor(x)
  levs <- levels(x)
  
  # Assign random 'weights' to each level that sum to the desired proportion
  probs <- stats::runif(length(levs))
  probs <- (probs / sum(probs)) * prop * length(levs)
  probs <- pmin(pmax(probs, 0), 1) # Ensure probabilities stay in [0, 1]
  
  # Map each individual observation to its level-based probability of being missing
  mask <- vapply(seq_along(x), function(i) {
    val <- x[i]
    if (is.na(val)) return(FALSE) 
    idx <- match(val, levs)
    stats::runif(1) < probs[idx]
  }, logical(1))
  
  return(mask)
}
