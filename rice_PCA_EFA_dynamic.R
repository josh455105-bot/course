# ================================================================
# Rice traits: data-driven MI -> PCA -> EFA
# No hard-coded variable names, missing columns, outliers, or # PCs
# User knows only: there should be 19 measured traits.
# Sources: supplied Factor Analysis lecture, supplied mi example.
# Requires R + packages "psych" and (only if NA exists) "mi".
# ================================================================

# [USER SETTINGS] You may change these settings after reviewing outputs.
datdir='/Users/ASUS/Downloads'
datfile=paste(datdir,'rice_wet.csv',sep='')      # If not found, RStudio opens file chooser
output_dir <- "rice_dynamic_results"
expected_traits <- 19L
trait_cols_override <- NULL     # Only if selection ambiguous (19 names)
mi_type_overrides <- character()  # Optional named values, e.g. c(trait_A="continuous")
manual_ncomp <- NULL           # NULL = use PCA parallel analysis suggestion
manual_nfactor <- NULL         # NULL = use EFA parallel analysis suggestion
mi_iterations <- 30L
mi_chains <- 4L
mi_rhat_limit <- 1.10
allow_poor_mi_convergence <- FALSE  # TRUE only after diagnosing convergence
mi_method_continuous <- "pmm"       # For missing numeric continuous variables
factor_extraction <- "minres"       # Alternative: "ml" (requires careful checks)
loading_threshold <- 0.50
parallel_repetitions <- 100L
set.seed(20261024)

if (!requireNamespace("psych", quietly = TRUE)) {
  stop('Please install psych first: install.packages("psych")')
}
if (!file.exists(datfile)) {
  if (interactive()) datfile <- file.choose()
  else stop("Input CSV not found: ", datfile)
}
if (!file.exists(datfile)) stop("Input file does not exist")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
write_out <- function(d, filename, row.names = FALSE) {
  utils::write.csv(d, file.path(output_dir, filename), row.names = row.names,
                   fileEncoding = "UTF-8", na = "")
}
section <- function(msg) cat("\n", strrep("=", 55), "\n", msg, "\n", sep = "")
`%||%` <- function(x, y) if (is.null(x)) y else x

# STEP 1. Read EVERYTHING as characters first (do not guess column types).
section("STEP 1 | Import, row count, column types")
missing_tokens <- c("", "NA", "N/A", "NaN", "NULL", ".", "?")
raw <- tryCatch(
  utils::read.csv(datfile, header = TRUE, check.names = FALSE,
                  colClasses = "character", stringsAsFactors = FALSE,
                  na.strings = missing_tokens, fileEncoding = "UTF-8-BOM"),
  error = function(e) stop("Cannot read CSV: ", conditionMessage(e),
                            ". Check file encoding / CSV delimiter."))
if (nrow(raw) < 3L || ncol(raw) < expected_traits) {
  stop("Not enough rows or columns. Rows: ", nrow(raw),
       "; columns: ", ncol(raw))
}
names(raw) <- make.unique(trimws(names(raw)))
raw[] <- lapply(raw, function(v) {
  v <- trimws(v)
  v[is.na(v) | v %in% missing_tokens] <- NA_character_
  v
})

# Numeric candidate = all PRESENT values parse as numeric. Never turn
# unexpected text into NA silently. All-NA columns are ambiguous.
numeric_candidate <- vapply(raw, function(v) {
  present <- v[!is.na(v)]
  length(present) > 0L && all(!is.na(suppressWarnings(as.numeric(present))))
}, logical(1))
profile <- data.frame(
  column = names(raw),
  numeric_candidate = unname(numeric_candidate),
  observed = vapply(raw, function(v) sum(!is.na(v)), integer(1)),
  missing = vapply(raw, function(v) sum(is.na(v)), integer(1)),
  unique_observed = vapply(raw, function(v) length(unique(stats::na.omit(v))), integer(1)),
  stringsAsFactors = FALSE
)
print(profile)
write_out(profile, "01_column_profile.csv")
cat("Rows:", nrow(raw), " | Columns:", ncol(raw), "\n")

if (is.null(trait_cols_override)) {
  selected <- names(raw)[numeric_candidate]
  if (length(selected) != expected_traits) {
    stop("Found ", length(selected), " numeric-like columns, not ",
         expected_traits, ". Read 01_column_profile.csv and supply",
         " trait_cols_override (19 names) before rerunning. No guessing.")
  }
} else {
  selected <- trait_cols_override
  if (length(selected) != expected_traits || anyDuplicated(selected) ||
      !all(selected %in% names(raw))) {
    stop("trait_cols_override must contain exactly 19 different existing columns")
  }
  if (!all(numeric_candidate[selected])) {
    stop("Manually selected trait columns contain nonnumeric text or",
         " entirely unobserved columns. Correct / document first.")
  }
}
metadata_cols <- setdiff(names(raw), selected)  # IDs are not assumed by name/position
metadata <- raw[, metadata_cols, drop = FALSE]
X <- as.data.frame(lapply(raw[, selected, drop = FALSE], as.numeric),
                   check.names = FALSE)
names(X) <- selected
if (any(!is.finite(as.matrix(X)) & !is.na(as.matrix(X)))) {
  stop("Infinite or non-finite observed values found. Check raw data.")
}
if (anyDuplicated(raw)) {
  warning("Duplicate entire raw rows exist. Please check; no rows removed.")
}
if (ncol(metadata) > 0L) {
  cat("Excluded metadata/non-numeric columns:",
      paste(names(metadata), collapse = ", "), "\n")
}
cat("Analysis variables:", paste(names(X), collapse = ", "), "\n")

# STEP 2. Report missingness and descriptive statistics. No deletion.
section("STEP 2 | Missingness, sample validity, possible outliers")
miss <- data.frame(variable = names(X),
                   n_missing = colSums(is.na(X)),
                   pct_missing = round(100 * colMeans(is.na(X)), 3),
                   n_unique = vapply(X, function(v) length(unique(stats::na.omit(v))), integer(1)))
print(miss)
write_out(miss, "02_missingness.csv")
cat("Rows with all traits missing:", sum(rowSums(!is.na(X)) == 0), "\n")
cat("Rows with at least one missing:", sum(!stats::complete.cases(X)), "\n")
if (any(rowSums(!is.na(X)) == 0)) {
  stop("One or more rows have ALL traits missing. Review these cases",
       " before MI. Rows were not removed automatically.")
}
if (any(colSums(!is.na(X)) < 5L)) {
  stop("One or more variables have <5 observed values. Review data.")
}

# Possible outliers, using 1.5 x IQR for FLAGGING only; never deleted.
flagged <- list()
for (v in names(X)) {
  z <- X[[v]]
  q <- stats::quantile(z, probs = c(0.25, 0.75), na.rm = TRUE)
  iq <- diff(q)
  if (is.finite(iq) && iq > 0) {
    lower <- q[1] - 1.5 * iq
    upper <- q[2] + 1.5 * iq
    idx <- which(!is.na(z) & (z < lower | z > upper))
    if (length(idx)) flagged[[v]] <- data.frame(
      row_number = idx, variable = v, observed_value = z[idx],
      lower_fence = lower, upper_fence = upper)
  }
}
outlier_flags <- if (length(flagged)) do.call(rbind, flagged) else
  data.frame(row_number = integer(), variable = character(),
             observed_value = numeric(), lower_fence = numeric(),
             upper_fence = numeric())
write_out(outlier_flags, "02_possible_outlier_flags_NO_DELETIONS.csv")
cat("Possible outlier records (not removed):", nrow(outlier_flags), "\n")
write_out(data.frame(row_number = seq_len(nrow(raw)), metadata, X,
                     check.names = FALSE), "02_raw_traits_and_metadata.csv")

if (any(vapply(X, function(v) stats::sd(v, na.rm = TRUE) == 0, logical(1)))) {
  stop("Constant observed trait exists. Inspect and decide what to do.")
}

# STEP 3. Only use MI if actually necessary; no assumed missing column.
section("STEP 3 | MI (if missing values exist)")
completed_sets <- NULL
if (!anyNA(X)) {
  cat("No missing trait values detected. MI skipped.\n")
  completed_sets <- list(X)
} else {
  if (!requireNamespace("mi", quietly = TRUE)) {
    stop('Missing values detected. Install mi: install.packages("mi")')
  }
  mdf <- mi::missing_data.frame(X)
  cat("MI's proposed data types and methods (REVIEW):\n")
  show(mdf)
  if (length(mi_type_overrides)) {
    if (is.null(names(mi_type_overrides)) ||
        !all(names(mi_type_overrides) %in% names(X))) {
      stop("mi_type_overrides must have names matching existing trait columns")
    }
    for (v in names(mi_type_overrides)) {
      mdf <- mi::change(mdf, y = v, what = "type",
                        to = unname(mi_type_overrides[[v]]))
    }
  }
  # Automatically apply PMM only to variables MI recognizes as continuous.
  # Do NOT assign ordinal/nominal types by column names or numeric ranges.
  missing_vars <- names(X)[colSums(is.na(X)) > 0]
  for (v in missing_vars) {
    var_obj <- mdf@variables[[v]]
    if (methods::is(var_obj, "continuous")) {
      mdf <- mi::change(mdf, y = v, what = "imputation_method",
                        to = mi_method_continuous)
    } else {
      cat("MI kept automatically detected noncontinuous method for:", v,
          " (class:", paste(class(var_obj), collapse = "/"), ")\n")
    }
  }
  cat("FINAL MI specifications - review before interpreting results:\n")
  show(mdf)
  capture.output(show(mdf), file = file.path(output_dir, "03_MI_variable_types.txt"))
  # If detected type is wrong, modify types using mi::change() and rerun.
  # MI assumptions (including missingness mechanism) cannot be tested by code alone.
  imp <- tryCatch(
    mi::mi(mdf, n.iter = mi_iterations, n.chains = mi_chains,
           max.minutes = Inf, parallel = FALSE),
    error = function(e) stop("MI failed; review predictors and model types: ",
                             conditionMessage(e)))
  rhat <- tryCatch(mi::Rhats(imp), error = function(e) {
    warning("R-hat diagnostic failed: ", conditionMessage(e)); NA_real_
  })
  print(rhat)
  write_out(data.frame(item = names(rhat) %||% seq_along(rhat),
                       Rhat = as.numeric(rhat)), "03_MI_Rhat.csv")
  finite_rhat <- as.numeric(rhat)[is.finite(as.numeric(rhat))]
  poor_convergence <- !length(finite_rhat) ||
    any(!is.finite(as.numeric(rhat))) ||
    any(finite_rhat > mi_rhat_limit)
  if (poor_convergence && !allow_poor_mi_convergence) {
    stop("MI convergence not verified (R-hat > ", mi_rhat_limit,
         " or unavailable). Extend/review MI, then rerun.")
  }
  completed_sets <- mi::complete(imp, m = mi_chains)
  if (!is.list(completed_sets)) completed_sets <- list(completed_sets)
  completed_sets <- lapply(completed_sets, function(d) {
    
    # 將原始欄名轉換成 R 合法欄名
    safe_names <- make.names(selected, unique = TRUE)
    
    # 確認 MI 完成資料中是否有這些欄位
    if (!all(safe_names %in% names(d))) {
      stop("Column names do not match after MI.")
    }
    
    # 只保留原本的 19 個性狀
    d <- d[, safe_names, drop = FALSE]
    
    # 恢復原始欄名
    names(d) <- selected
    
    # 轉換為數值
    d[] <- lapply(d, function(v) as.numeric(as.character(v)))
    
    # 確認沒有剩餘缺失值
    if (anyNA(d)) {
      stop("NA remains after imputation. Review MI output.")
    }
    
    d
  })
  
}

# No averaging of imputed cells; keep completed datasets separate.
for (i in seq_along(completed_sets)) {
  d <- completed_sets[[i]]
  write_out(data.frame(row_number = seq_len(nrow(d)), metadata, d,
                       check.names = FALSE),
            sprintf("03_completed_dataset_%02d.csv", i))
}

# STEP 4. Data diagnostics BEFORE PCA / EFA, using completed set 1.
section("STEP 4 | Descriptive stats, correlations, Bartlett, KMO")
Y <- completed_sets[[1]]
n <- nrow(Y)
p <- ncol(Y)
desc <- psych::describe(Y)
write_out(data.frame(variable = rownames(desc), desc, row.names = NULL),
          "04_descriptive_statistics.csv")
R <- stats::cor(Y)
write_out(data.frame(variable = rownames(R), R, check.names = FALSE),
          "04_correlation_matrix.csv")
possible_collinearity <- which(abs(R) >= 0.95 & upper.tri(R), arr.ind = TRUE)
if (nrow(possible_collinearity)) {
  high_cor <- data.frame(var1 = rownames(R)[possible_collinearity[, 1]],
                         var2 = colnames(R)[possible_collinearity[, 2]],
                         correlation = R[possible_collinearity])
} else {
  high_cor <- data.frame(var1 = character(), var2 = character(),
                         correlation = numeric())
}
write_out(high_cor, "04_very_high_correlations_review.csv")
if (nrow(high_cor)) warning("Some correlations >= |0.95|; review redundancy")
if (min(eigen(R, symmetric = TRUE, only.values = TRUE)$values) <= 1e-8) {
  stop("Correlation matrix singular/nearly singular. Check redundant or",
       " derived traits before PCA/EFA; nothing deleted automatically.")
}
print(psych::KMO(R))
print(psych::cortest.bartlett(R, n = n))
kmo <- psych::KMO(R)
bart <- psych::cortest.bartlett(R, n = n)
write_out(data.frame(variable = names(kmo$MSAi), MSA = unname(kmo$MSAi)),
          "04_MSA_by_variable.csv")
write_out(data.frame(KMO_overall = unname(kmo$MSA),
                     Bartlett_chisq = unname(bart$chisq),
                     Bartlett_df = unname(bart$df),
                     Bartlett_p = unname(bart$p.value)),
          "04_KMO_Bartlett_summary.csv")
if (kmo$MSA < 0.6) warning("KMO < 0.6. EFA factors may be difficult to justify")

# STEP 5. Discover dimension counts, do NOT fix five components/factors.
section("STEP 5 | Eigenvalues / scree / parallel analysis")
pca <- stats::prcomp(Y, center = TRUE, scale. = TRUE)
eig <- pca$sdev^2
ev <- data.frame(PC = paste0("PC", seq_len(p)), eigenvalue = eig,
                 explained_pct = 100 * eig / sum(eig),
                 cumulative_pct = 100 * cumsum(eig) / sum(eig))
print(ev)
write_out(ev, "05_PCA_eigenvalues.csv")
png(file.path(output_dir, "05_scree_plot.png"), width = 1000, height = 700)
plot(seq_len(p), eig, type = "b", pch = 19, xlab = "Component",
     ylab = "Eigenvalue", main = "Unrotated PCA scree plot")
abline(h = 1, col = "red", lty = 2)
dev.off()

grDevices::png(file.path(output_dir, "05_parallel_analysis.png"),
               width = 1100, height = 800)
pa <- tryCatch(
  psych::fa.parallel(Y, fm = factor_extraction, fa = "both",
                     n.iter = parallel_repetitions, plot = TRUE),
  error = function(e) {
    warning("Parallel analysis failed: ", conditionMessage(e)); NULL
  })
grDevices::dev.off()
if (is.null(pa)) {
  if (is.null(manual_ncomp) || is.null(manual_nfactor)) {
    stop("Parallel analysis failed. Inspect the data and, if scientifically",
         " justified, supply BOTH manual_ncomp and manual_nfactor.")
  }
  pa <- list(ncomp = NA_integer_, nfact = NA_integer_)
}
cat("PA suggested PCs:", pa$ncomp, "; factors:", pa$nfact, "\n")
write_out(data.frame(method = c("PCA_parallel", "EFA_parallel",
                                "PCA_Kaiser", "PCA_to_80pct"),
                     suggested_count = c(pa$ncomp, pa$nfact,
                                         sum(eig > 1),
                                         which(ev$cumulative_pct >= 80)[1])),
          "05_dimension_decisions.csv")
k_pc <- if (is.null(manual_ncomp)) pa$ncomp else manual_ncomp
k_f <- if (is.null(manual_nfactor)) pa$nfact else manual_nfactor
valid_k <- function(k) length(k) == 1L && is.finite(k) && k == as.integer(k) &&
  k >= 1L && k < p
if (!valid_k(k_pc) || !valid_k(k_f)) {
  stop("Suggested component/factor count invalid (e.g. 0, NA).",
       " Review diagnostics; set manual_ncomp/manual_nfactor if justified.")
}
k_pc <- as.integer(k_pc)
k_f <- as.integer(k_f)

# STEP 6. Raw PCA (weights are NOT factor loadings).
section("STEP 6 | Raw PC weights / loadings / contributions")
W <- pca$rotation[, seq_len(k_pc), drop = FALSE]
for (j in seq_len(k_pc)) {
  if (W[which.max(abs(W[, j])), j] < 0) W[, j] <- -W[, j]
}
L_raw <- sweep(W, 2, pca$sdev[seq_len(k_pc)], FUN = "*")
contr <- 100 * W^2
write_out(data.frame(variable = rownames(W), W, check.names = FALSE),
          "06_PCA_raw_weights.csv")
write_out(data.frame(variable = rownames(L_raw), L_raw, check.names = FALSE),
          "06_PCA_raw_loadings.csv")
write_out(data.frame(variable = rownames(contr), contr, check.names = FALSE),
          "06_PCA_contribution_percent.csv")
ranked <- do.call(rbind, lapply(seq_len(k_pc), function(j) {
  ord <- order(abs(W[, j]), decreasing = TRUE)
  data.frame(PC = paste0("PC", j), rank = seq_along(ord),
             variable = rownames(W)[ord], weight = W[ord, j],
             loading = L_raw[ord, j], contribution_pct = contr[ord, j])
}))
write_out(ranked, "06_PCA_weights_ranked_by_PC.csv")
print(subset(ranked, rank <= 5))

# STEP 7. Rotated PCA, chosen for interpretability (not identical to raw PCs).
section("STEP 7 | Rotated PCA components")
pca_rot <- psych::principal(Y, nfactors = k_pc,
                            rotate = if (k_pc > 1L) "varimax" else "none")
L_rc <- as.matrix(unclass(pca_rot$loadings))
write_out(data.frame(variable = rownames(L_rc), L_rc, check.names = FALSE),
          "07_PCA_rotated_loadings.csv")
print(pca_rot$loadings, cutoff = 0.30)

# STEP 8. EFA, WITHOUT forcing factor count to match PCA count.
section("STEP 8 | EFA factor loadings, communality and uniqueness")
efa <- tryCatch(
  psych::fa(Y, nfactors = k_f, rotate = if (k_f > 1L) "varimax" else "none",
            fm = factor_extraction, scores = "regression"),
  error = function(e) stop("EFA failed: ", conditionMessage(e),
                           ". Check KMO, redundancy, and factor count."))
print(efa$loadings, cutoff = 0.30)
L_e <- as.matrix(unclass(efa$loadings))
write_out(data.frame(variable = rownames(L_e), L_e, check.names = FALSE),
          "08_EFA_loadings.csv")
property <- data.frame(variable = rownames(L_e),
                       communality = as.numeric(efa$communality),
                       uniqueness = as.numeric(efa$uniquenesses),
                       max_abs_loading = apply(abs(L_e), 1, max),
                       n_factors_ge_threshold = rowSums(abs(L_e) >= loading_threshold))
write_out(property, "08_EFA_communality_uniqueness.csv")
count_table <- data.frame(factor = colnames(L_e),
                          n_loading_abs_ge_threshold = colSums(abs(L_e) >= loading_threshold),
                          SS_loadings = colSums(L_e^2),
                          proportion_of_total_variance = colSums(L_e^2) / p)
write_out(count_table, "08_EFA_factor_summary.csv")
print(count_table)
if (any(efa$uniquenesses < 0.01)) warning("Very low uniqueness; possible Heywood/redundant traits")
if (any(count_table$n_loading_abs_ge_threshold < 2)) {
  warning("Factor supported by <2 strong-loading variables; interpret cautiously")
}
if (any(property$n_factors_ge_threshold > 1)) {
  warning("Cross-loading variables found; inspect 08_EFA_communality_uniqueness.csv")
}
if (!is.null(efa$Vaccounted)) {
  write_out(data.frame(measure = rownames(efa$Vaccounted),
                       efa$Vaccounted, check.names = FALSE),
            "08_EFA_variance_accounted.csv")
}

# STEP 9. PCA/EFA pattern comparison, including unequal dimension cases.
section("STEP 9 | PCA vs EFA congruence")
congruence <- psych::factor.congruence(L_rc, L_e)
print(congruence)
write_out(data.frame(rotated_PC = rownames(congruence), congruence,
                     check.names = FALSE), "09_congruence_matrix.csv")
cat("NOTE: PC/Factor order and signs can differ. Compare the FULL matrix.\n")

# STEP 10. Figures and across-imputation sensitivity.
section("STEP 10 | Plots and multiple-imputation sensitivity")
plot_top <- function(m, k, file) {
  num_show <- k
  png(file, width = 1300, height = 350 * ceiling(num_show / 2))
  op <- par(mfrow = c(ceiling(num_show / 2), 2), mar = c(4, 10, 4, 2))
  on.exit({par(op); dev.off()}, add = TRUE)
  for (j in seq_len(num_show)) {
    sel <- order(abs(m[, j]), decreasing = TRUE)[seq_len(min(7L, nrow(m)))]
    vals <- rev(m[sel, j])
    barplot(vals, horiz = TRUE, las = 1,
            names.arg = rev(rownames(m)[sel]),
            col = ifelse(vals >= 0, "steelblue", "tomato"),
            main = colnames(m)[j], xlab = "Signed weight / loading")
    abline(v = 0, col = "gray50")
  }
}
plot_top(W, k_pc, file.path(output_dir, "10_PCA_top_raw_weights.png"))
plot_top(L_e, k_f, file.path(output_dir, "10_EFA_top_factor_loadings.png"))

# Factor scores yield a FREQUENCY DISTRIBUTION of samples; not factor importance.
if (!is.null(efa$scores)) {
  scores <- as.data.frame(efa$scores)
  write_out(data.frame(row_number = seq_len(nrow(scores)), scores),
            "10_EFA_factor_scores.csv")
  png(file.path(output_dir, "10_EFA_factor_score_histograms.png"),
      width = 1300, height = 350 * ceiling(ncol(scores) / 2))
  op <- par(mfrow = c(ceiling(ncol(scores) / 2), 2), mar = c(4, 4, 4, 2))
  for (j in seq_len(ncol(scores))) {
    hist(scores[[j]], breaks = "FD", main = colnames(scores)[j],
         xlab = "Estimated factor score", col = "gray80", border = "white")
    abline(v = 0, col = "red", lty = 2)
  }
  par(op)
  dev.off()
}

# Do not pool rotated loadings by simply taking arithmetic averages; the
# order/sign of factors can change across imputed datasets. Instead, report
# PCA eigenvalue stability and agreement of factor patterns with set 1.
if (length(completed_sets) > 1L) {
  pca_robust <- do.call(rbind, lapply(seq_along(completed_sets), function(i) {
    cur <- stats::prcomp(completed_sets[[i]], center = TRUE, scale. = TRUE)
    e <- cur$sdev^2
    data.frame(imputation = i, PC = seq_along(e), eigenvalue = e,
               explained_pct = 100 * e / sum(e))
  }))
  write_out(pca_robust, "10_MI_PCA_eigenvalue_sensitivity.csv")
  # Rotated factors in each set should be compared with the reference
  # using the full congruence matrix, not matched by their raw names.
  efa_robust <- list()
  for (i in seq_along(completed_sets)) {
    ei <- tryCatch(
      psych::fa(completed_sets[[i]], nfactors = k_f,
                rotate = if (k_f > 1) "varimax" else "none",
                fm = factor_extraction),
      error = function(e) NULL)
    if (!is.null(ei)) {
      ci <- psych::factor.congruence(L_e,
                                     as.matrix(unclass(ei$loadings)))
      temp <- as.data.frame(as.table(ci))
      names(temp) <- c("reference_factor", "comparison_factor", "congruence")
      temp$imputation <- i
      efa_robust[[length(efa_robust) + 1L]] <- temp
    } else warning("EFA failed on imputed set ", i, "; inspect stability")
  }
  if (length(efa_robust)) {
    write_out(do.call(rbind, efa_robust), "10_MI_EFA_congruence_sensitivity.csv")
  }
}

section("FINISHED")
cat("All outputs stored in: ", normalizePath(output_dir), "\n", sep = "")
cat("No traits or observations automatically removed.\n")
cat("Dimension counts selected from data, not fixed at five.\n")
cat("Multiple-imputation results are exploratory; these EFA loadings are\n")
cat("NOT a formal Rubin-rules pooled EFA inference.\n")
cat("Please interpret all factors using agronomic knowledge before naming.\n")

