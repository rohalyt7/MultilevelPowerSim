# **MultilevelPowerSim: Simulation-Based Power Analysis**
**Author:** Thomas Rohaly  
**Role:** Computational Research Psychologist  
**Status:** Active Development

## **Overview**
Statistical power in multilevel modeling (MLM) is notoriously difficult to calculate using analytic formulas due to the complex interplay between sample size ($N$), cluster size ($n$), Intraclass Correlation (ICC), and random effect structures. 

This repository provides a custom, modular framework to conduct **Simulation-Based Power Analysis**. By simulating thousands of datasets with known "ground truth" parameters, researchers can empirically determine the power, bias, and convergence rates of their models before data collection begins.

## **Core Philosophy: Simulation-Driven Rigor**
Consistent with a modular ETL approach, this project emphasizes Data Simulation as a prerequisite for robust analysis. By building custom research frameworks, we can test models against known ground truths before applying them to "noisy" real-world physiological and behavioral data streams.

## **Current Modules**
The framework is organized into a logical dependency flow to ensure reproducibility:

### **1. Simulation Core (`R/01_sim_core.R`)**
The foundation for generating hierarchical data.
- **Functionality:** Supports one-shot generation of clustered data or modular expansion of fixed-effects datasets into multilevel structures.
- **Key Features:** Custom random intercept and random slope injection using `MASS::mvrnorm`.

### **2. Data Modifiers (`R/02_modifiers.R`)**
Tools to inject real-world complexity into synthetic data.
- **Capabilities:** Adding level-2 covariates, generating longitudinal time-points, forcing predictor correlations, and inducing heteroskedasticity.
- **Statistical Rigor:** Includes advanced centering options (Person-Mean Centering) to isolate within-person vs. between-person effects.

### **3. Missing Data Engine (`R/03_missing_data.R`)**
A robust suite for testing model sensitivity to data loss—a critical concern in human subjects research.
- **Mechanisms:** Supports **MCAR** (Random), **MAR** (Conditional on covariates), and **MNAR** (Conditional on the missing value itself).

### **4. Power Engine & Visualization (`R/04_power_engine.R`)**
The "Brain" of the framework that iterates through simulations and handles model fitting via `lme4`.
- **Diagnostics:** Automatically tracks and reports singularity (over-parameterization) and convergence warnings.
- **Visualization:** Integrated plotting for Power Curves and Power Heatmaps across varying Effect Sizes and ICC levels.

---

## **Technical Notes & Statistical Caveats**
To ensure the integrity of the results, users should be aware of the following design considerations implemented in this framework:

* **Convergence & Singularity:** Mixed-effects models often result in "Singular Fits" (e.g., when random effect variances are estimated at zero). This engine explicitly flags these instances rather than ignoring them, allowing researchers to assess model stability.
* **Stochastic vs. Deterministic MAR:** The `add_missing_data` module utilizes a quantile-based threshold for MAR/MNAR. This provides a clear, interpretable "Top-N%" missingness structure.
* **The ICC/Slope Variance Trade-off:** When adding random slopes, the total variance in the model increases. The ICC calculation provided is based on the intercept-only portion; researchers should be mindful of how random slopes impact total variance.
* **Handling Categorical Slopes:** To prevent mathematical errors during simulation, categorical predictors used as random slopes are automatically cast to numeric vectors. 

---

## **Roadmap & Future Additions**
- **Unit Testing Suite:** Development of formal test scripts for each module to verify that data generation and manipulation functions remain robust across edge cases.
- **GEE Integration:** Expanding the power engine to support Generalized Estimating Equations (GEE) as an alternative for clustered data.
- **Bayesian Frameworks:** Transitioning the fit functions to support `brms` for robust hierarchical modeling.
- **Example Usage Library:** Generation of specialized scripts demonstrating power analyses for complex longitudinal designs and physiological signal "dropouts."

## **Getting Started**
1. Clone the repository.
2. Open `MultilevelPowerSim.Rproj`.
3. Load all functions: `lapply(list.files("R", full.names = TRUE), source)`
4. (Coming Soon) See `scripts/example_usage.R` for a full pipeline walkthrough.

## **License**
This project is licensed under the MIT License - see the LICENSE file for details.
