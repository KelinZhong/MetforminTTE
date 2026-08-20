# ---------------------------------------------------------------------------
# The pipeline: validate, then run.
#
# THREE ORDERINGS THAT MUST NOT CHANGE.
#
# 1. Matching (.tte_match_impl): caliper prefilter -> fit the propensity score
#    ON THE PREFILTERED SAMPLE -> greedy 1:1 match under the joint caliper.
#    Fitting on the full cohort changes every score and every pair.
#
# 2. Analysis time (tte_analysis_time): set time zero -> apply the strategy ->
#    derive outcomes. Always that order. Derivation reads the already-rewritten
#    incidence and date columns; deriving first and censoring second produces a
#    different analysis with no error and no warning.
#
# 3. Time zero: both members of a pair take the EXPOSED member\'s date. This is
#    why derivation must follow matching rather than being precomputed.
# ---------------------------------------------------------------------------

utils::globalVariables("subclass")

#' Validate inputs before running anything
#'
#' Checks every distinct source (not just the first — files built at different
#' times can diverge in columns or coding) and prints a frequency table for
#' every prevalence, status and incidence column.
#'
#' The tables are the point. A prevalence flag levelled `"1"`/`"2"` rather than
#' 0/1 will be compared against 1 downstream and silently select the
#' *unaffected* group in every exclusion — an error that produces plausible
#' output and no warning. Printing the coding on every run makes it visible
#' instead of latent; feed the answer back through
#' `tte_params(prev_affected_level = )`.
#'
#' @param x A [tte_scenarios()] object, a character vector of paths, or a
#'   data.frame.
#' @param outcomes An outcome registry from [tte_outcomes()].
#' @param params A [tte_params()] object.
#' @param rx Optional dispensing records to validate as well: reports
#'   unparseable dates, duplicates, and how many exposed participants have no
#'   records at all.
#' @param eid_col,date_col Column names in the dispensing data.
#' @param verbose Print the tables. Default `TRUE`.
#'
#' @return Invisibly, a list with the cohort summary and, if `rx` was given,
#'   the dispensing summary.
#' @examples
#' set.seed(1); n <- 200
#' cohort <- data.frame(
#'   pid = sprintf("P%03d", seq_len(n)),
#'   exposed = rep(c(1L, 0L), c(60, 140)),
#'   age = round(stats::rnorm(n, 60, 8)), biom = stats::rnorm(n, 6, 0.6),
#'   base_date = as.Date("2010-01-01"),
#'   dth = stats::rbinom(n, 1, 0.2), dth_date = as.Date("2015-06-01"),
#'   oa_prev = stats::rbinom(n, 1, 0.1), oa_status = 0L,
#'   oa_incid = stats::rbinom(n, 1, 0.3), oa_date = as.Date("2014-03-01"))
#' outcomes <- tte_outcomes("oa", "Osteoarthritis", "oa_prev", "oa_status",
#'                          "oa_incid", "oa_date", model = "cox")
#' params <- tte_params("pid", "exposed", exposed ~ age + biom, "base_date",
#'                      death_status_col = "dth", death_date_col = "dth_date",
#'                      match_calipers = c(biom = 0.5), ps_caliper = 0.2,
#'                      followup_end = as.Date("2019-12-31"))
#'
#' info <- tte_validate(cohort, outcomes, params)
#' info$cohort[[1]]$tables
#' @export
tte_validate <- function(x, outcomes, params, rx = NULL,
                         eid_col = "eid", date_col = "issue_date",
                         verbose = TRUE) {
  res <- list(cohort = .validate_cohort(x, outcomes, params, verbose))
  if (!is.null(rx)) {
    src <- if (inherits(x, "tte_scenarios")) x$source[1] else
      if (is.list(x) && !is.data.frame(x)) x[[1]] else x
    res$rx <- .validate_rx(rx, src, params, eid_col, date_col, verbose)
  }
  invisible(res)
}
#' @noRd
.validate_cohort <- function(x, outcomes, params, verbose = TRUE) {

  sources <- if (inherits(x, "tte_scenarios")) unique(x$source) else x
  if (is.data.frame(sources)) sources <- list(sources)

  res <- list()
  for (i in seq_along(sources)) {
    src <- if (is.list(sources)) sources[[i]] else sources[i]
    nm <- if (is.character(src)) src else paste0("<data.frame ", i, ">")
    if (verbose) message("\n=== ", nm, " ===")

    data <- tte_read(src, id_col = params$id_col)
    missing_cols <- require_cols(data, outcomes, params, error = FALSE)

    elem <- outcomes[outcomes$model != "derived", ]
    cols <- unique(stats::na.omit(c(elem$prev, elem$status, elem$incid)))
    cols <- intersect(cols, names(data))

    tabs <- lapply(cols, function(cl) table(data[[cl]], useNA = "ifany"))
    names(tabs) <- cols

    if (verbose) {
      if (length(missing_cols)) {
        message("MISSING COLUMNS: ", paste(missing_cols, collapse = ", "))
      } else {
        message("All required columns present (", ncol(data), " columns, ",
                nrow(data), " rows).")
      }
      message("Exposure: ",
              paste(utils::capture.output(
                table(data[[params$exposure_col]], useNA = "ifany")),
                collapse = " | "))
      message("--- prevalence / status / incidence coding ---")
      for (cl in cols) {
        message("  ", cl, ": ",
                paste(names(tabs[[cl]]), unname(tabs[[cl]]),
                      sep = "=", collapse = "  "))
      }
    }

    res[[nm]] <- list(missing = missing_cols, tables = tabs,
                      nrow = nrow(data), ncol = ncol(data))
  }

  invisible(res)
}
#' @noRd
.validate_rx <- function(rx, cohort = NULL, params,
                         eid_col = "eid", date_col = "issue_date",
                         verbose = TRUE) {

  d <- tte_read(rx)
  miss <- setdiff(c(eid_col, date_col), names(d))
  if (length(miss)) {
    stop("Dispensing file missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }

  parsed <- .as_date(d[[date_col]])
  n_unparsed <- sum(is.na(parsed) & !is.na(d[[date_col]]))
  examples <- utils::head(unique(d[[date_col]][is.na(parsed) &
                                                 !is.na(d[[date_col]])]), 5)

  dup <- sum(duplicated(d[, c(eid_col, date_col), with = FALSE]))

  n_zero_rx <- NA_integer_
  if (!is.null(cohort)) {
    ch <- tte_read(cohort, id_col = params$id_col)
    exposed <- ch[[params$id_col]][ch[[params$exposure_col]] == 1]
    n_zero_rx <- sum(!exposed %in% unique(d[[eid_col]]))
  }

  if (verbose) {
    message("Dispensing records: ", nrow(d), " rows, ",
            length(unique(d[[eid_col]])), " unique participants.")
    message("Unparseable dates: ", n_unparsed,
            if (n_unparsed) paste0(" (e.g. ", paste(examples, collapse = ", "),
                                   ")") else "")
    message("Exact duplicate rows: ", dup)
    if (!is.na(n_zero_rx)) {
      message("Exposed participants with zero dispensing records: ", n_zero_rx,
              " -- see `missing_rx` in tte_strategy().")
    }
  }

  invisible(list(n_rows = nrow(d), n_unparsed = n_unparsed,
                 unparsed_examples = examples, n_duplicates = dup,
                 n_zero_rx = n_zero_rx))
}
#' @noRd
.load_cohort <- function(source, outcomes, params) {
  data <- tte_read(source, id_col = params$id_col,
                   reader = params$reader)
  data <- params$recode(data)
  data <- data.table::as.data.table(data)
  require_cols(data, outcomes, params)
  if (!is.null(params$prev_affected_level)) {
    data <- normalize_prev(data, outcomes, params$prev_affected_level)
  }
  data
}
#' @noRd
normalize_prev <- function(data, outcomes, affected_level, cols = NULL) {

  if (missing(affected_level) || length(affected_level) != 1L) {
    stop("`affected_level` is required: read it off the tables printed by ",
         "tte_validate(). Do not guess.", call. = FALSE)
  }
  data <- data.table::as.data.table(data)
  if (is.null(cols)) {
    cols <- unique(stats::na.omit(outcomes$prev))
  }
  cols <- intersect(cols, names(data))

  for (cl in cols) {
    v <- as.character(data[[cl]])
    lv <- setdiff(unique(v), NA)
    if (!as.character(affected_level) %in% lv) {
      warning("Column \"", cl, "\" has no level \"", affected_level,
              "\" (levels: ", paste(lv, collapse = ", "),
              "); result will be all zero.", call. = FALSE)
    }
    data[[cl]] <- as.integer(v == as.character(affected_level))
  }
  data
}
#' @noRd
tte_read <- function(path, id_col = NULL, reader = NULL, ...) {

  if (is.data.frame(path)) {
    out <- data.table::as.data.table(path)
    return(.sort_by_id(out, id_col))
  }
  if (!is.character(path) || length(path) != 1L) {
    stop("`path` must be a single file path or a data.frame.", call. = FALSE)
  }
  if (!file.exists(path)) {
    stop("File not found: ", path, call. = FALSE)
  }

  if (!is.null(reader)) {
    out <- reader(path, ...)
    return(.sort_by_id(data.table::as.data.table(out), id_col))
  }

  ext <- tolower(tools::file_ext(path))
  out <- switch(
    ext,
    csv = data.table::fread(path, ...),
    tsv = data.table::fread(path, ...),
    txt = data.table::fread(path, ...),
    rds = readRDS(path),
    parquet = stop("Reading .parquet directly is not supported. Pass the ",
                   "reader explicitly:\n  tte_read(path, reader = ",
                   "arrow::read_parquet)\nor convert the file to .csv or ",
                   ".rds first.", call. = FALSE),
    stop("Unsupported extension \".", ext, "\". Pass `reader` to override.",
         call. = FALSE)
  )

  .sort_by_id(data.table::as.data.table(out), id_col)
}
#' @noRd
.sort_by_id <- function(dt, id_col) {
  if (is.null(id_col) || !id_col %in% names(dt)) return(dt)
  dt[order(dt[[id_col]]), ]
}
#' @noRd
tte_recode_ukb <- function(data, education_col = "education",
                           ethnicity_col = "ethnicity") {

  data <- data.table::as.data.table(data)

  if (education_col %in% names(data)) {
    e <- as.character(data[[education_col]])
    data[[education_col]] <- factor(
      data.table::fcase(
        e %in% c("1", "College or University degree"), "Degree",
        e %in% c("2", "3", "4", "5", "6"), "Post-secondary or vocational",
        e %in% c("-7", "None of the above"), "None of the above",
        default = NA_character_),
      levels = c("None of the above", "Post-secondary or vocational", "Degree"))
  }

  if (ethnicity_col %in% names(data)) {
    x <- as.character(data[[ethnicity_col]])
    data[[ethnicity_col]] <- factor(
      data.table::fcase(
        substr(x, 1L, 1L) == "1" | grepl("^White", x), "White",
        substr(x, 1L, 1L) == "2" | grepl("Mixed", x), "Mixed",
        substr(x, 1L, 1L) == "3" | grepl("Asian", x), "Asian",
        substr(x, 1L, 1L) == "4" | grepl("Black", x), "Black",
        default = "Other"),
      levels = c("White", "Mixed", "Asian", "Black", "Other"))
  }

  data
}
#' @noRd
tte_match <- function(cohort, params, cache_dir = NULL, source = NULL,
                      force = FALSE) {

  if (is.null(source)) source <- attr(cohort, "source")
  compute <- function() .tte_match_impl(cohort, params)

  ## Caching keys on the source file's identity and content signature, so an
  ## in-memory data.frame cannot be cached; fall through rather than failing
  ## inside normalizePath().
  if (is.null(cache_dir)) return(compute())
  if (!is.character(source) || length(source) != 1L) {
    message("tte_match(): caching skipped -- `source` is not a file path.")
    return(compute())
  }
  tte_cache(source, params, cache_dir, compute, force = force)
}
#' @noRd
.tte_match_impl <- function(cohort, params) {

  data <- data.table::as.data.table(cohort)
  exp_col <- params$exposure_col
  id_col <- params$id_col

  is_user <- !is.na(data[[exp_col]]) & data[[exp_col]] == 1
  u_idx <- which(is_user)
  c_idx <- which(!is_user)
  if (!length(u_idx) || !length(c_idx)) {
    stop("Need both exposed and unexposed participants; found ",
         length(u_idx), " and ", length(c_idx), ".", call. = FALSE)
  }

  ## --- step 1: caliper prefilter ------------------------------------------
  if (as.numeric(length(u_idx)) * length(c_idx) > 2e8) {
    message("Eligibility matrix is ", length(u_idx), " x ", length(c_idx),
            " -- this is memory-intensive.")
  }
  elig <- matrix(TRUE, nrow = length(u_idx), ncol = length(c_idx))
  for (cl in names(params$match_calipers)) {
    cal <- params$match_calipers[[cl]]
    uv <- as.numeric(data[[cl]])[u_idx]
    cv <- as.numeric(data[[cl]])[c_idx]
    d <- abs(outer(uv, cv, "-"))
    d[is.na(d)] <- Inf
    elig <- elig & (d <= cal)
  }

  keep_u <- rowSums(elig) > 0
  keep_c <- colSums(elig) > 0
  u_idx <- u_idx[keep_u]
  c_idx <- c_idx[keep_c]
  elig <- elig[keep_u, keep_c, drop = FALSE]
  if (!length(u_idx)) {
    stop("No exposed participant has an eligible control after the caliper ",
         "prefilter. Check `match_calipers`.", call. = FALSE)
  }

  pre <- data[sort(c(u_idx, c_idx)), ]

  ## --- step 2: propensity score, fitted on the prefiltered sample ---------
  fit <- stats::glm(params$ps_formula, data = pre,
                    family = stats::quasibinomial())
  pre_ps <- stats::predict(fit, newdata = pre, type = "response")
  ps_by_id <- stats::setNames(as.numeric(pre_ps), as.character(pre[[id_col]]))

  ps_u <- ps_by_id[as.character(data[[id_col]][u_idx])]
  ps_c <- ps_by_id[as.character(data[[id_col]][c_idx])]

  ## --- step 3: greedy 1:1 matching ----------------------------------------
  m <- greedy_match(elig = elig, ps_u = ps_u, ps_c = ps_c,
                    ps_caliper = params$ps_caliper)

  matched_u <- which(!is.na(m$control))
  if (!length(matched_u)) {
    stop("No pairs formed. Check `ps_caliper` and `match_calipers`.",
         call. = FALSE)
  }

  user_rows <- u_idx[matched_u]
  ctrl_rows <- c_idx[m$control[matched_u]]
  subclass <- seq_along(matched_u)

  out <- rbind(data[user_rows, ], data[ctrl_rows, ])
  out$subclass <- c(subclass, subclass)
  out$ps <- c(ps_u[matched_u], ps_c[m$control[matched_u]])
  out <- out[order(out$subclass, -out[[exp_col]]), ]

  structure(
    list(data = out,
         match_matrix = elig,
         match_num = rowSums(elig),
         user_ids = data[[id_col]][user_rows],
         control_ids = data[[id_col]][ctrl_rows],
         ps = ps_by_id,
         ps_model = fit,
         n_prefiltered = nrow(pre),
         cache_key = NA_character_,
         params_used = .matching_params(params)),
    class = "tte_matched_set")
}
#' @noRd
greedy_match <- function(elig, ps_u, ps_c, ps_caliper) {

  n_u <- nrow(elig)
  used <- rep(FALSE, ncol(elig))
  control <- rep(NA_integer_, n_u)
  distance <- rep(NA_real_, n_u)

  order_u <- order(rowSums(elig), method = "radix")

  for (i in order_u) {
    cand <- which(elig[i, ] & !used)
    if (!length(cand)) next
    d <- abs(ps_c[cand] - ps_u[i])
    d[is.na(d)] <- Inf
    ok <- d <= ps_caliper
    if (!any(ok)) next
    cand <- cand[ok]; d <- d[ok]
    pick <- cand[which.min(d)]
    control[i] <- pick
    distance[i] <- min(d)
    used[pick] <- TRUE
  }

  list(control = control, distance = distance)
}
#' @noRd
tte_analysis_time <- function(matched, outcomes,
                              strategy = tte_strategy("itt"), params) {

  if (!inherits(matched, "tte_matched_set")) {
    stop("`matched` must be a tte_matched_set from tte_match().", call. = FALSE)
  }
  if (!inherits(strategy, "tte_strategy")) {
    stop("`strategy` must be a tte_strategy object.", call. = FALSE)
  }

  surv <- set_time_zero(matched, params)
  attr(surv, "outcomes") <- outcomes

  surv <- strategy$process(surv, matched, params)

  ## processor contract: derived columns must not exist yet
  derived_now <- grep("^time_to_|_incid_cr$", names(surv), value = TRUE)
  if (length(derived_now)) {
    stop("Strategy \"", strategy$label, "\" returned derived column(s): ",
         paste(derived_now, collapse = ", "),
         ". A strategy processor may rewrite status / incidence / date only; ",
         "derivation happens after it.", call. = FALSE)
  }

  surv <- derive_outcomes(surv, outcomes, params)

  attr(surv, "strategy_label") <- strategy$label
  attr(surv, "followup_end") <- params$followup_end
  attr(surv, "time_zero_rule") <-
    paste0("exposed member's `", params$time_zero_col,
           "`, applied to both members of the pair")
  class(surv) <- unique(c("tte_surv", class(surv)))
  surv
}
#' @noRd
set_time_zero <- function(matched, params) {

  d <- data.table::as.data.table(matched$data)
  exp_col <- params$exposure_col
  exposed <- !is.na(d[[exp_col]]) & d[[exp_col]] == 1

  t0_user <- stats::setNames(
    as.numeric(.as_date(d[[params$time_zero_col]])[exposed]),
    as.character(d[["subclass"]][exposed]))

  missing_sc <- setdiff(as.character(d[["subclass"]]), names(t0_user))
  if (length(missing_sc)) {
    stop(length(missing_sc), " matched set(s) have no exposed member; ",
         "cannot assign time zero.", call. = FALSE)
  }

  d[["time_zero"]] <- as.Date(unname(t0_user[as.character(d[["subclass"]])]),
                              origin = "1970-01-01")
  d
}
#' @noRd
derive_outcomes <- function(surv, outcomes, params) {

  ## A strategy processor may write a `censor_date` column; follow-up is capped
  ## at it. This is the contract that keeps per-protocol person-time from
  ## exceeding intention-to-treat person-time.
  cens <- if ("censor_date" %in% names(surv)) surv[["censor_date"]] else NULL

  elem <- outcomes[outcomes$model != "derived", ]
  for (i in seq_len(nrow(elem))) {
    use_cr <- !identical(elem$model[i], "cox")
    surv <- derive_outcome(
      surv,
      key = elem$key[i],
      incid = elem$incid[i],
      date = elem$date[i],
      year_days = params$year_days,
      followup_end = params$followup_end,
      death_status = if (use_cr) params$death_status_col else NULL,
      death_date = if (use_cr) params$death_date_col else NULL,
      censor_date = cens)
  }

  der <- outcomes[outcomes$model == "derived", ]
  for (i in seq_len(nrow(der))) {
    surv <- derive_composite_endpoint(surv, der[i, ], outcomes, params)
  }

  surv
}
#' @noRd
tte_analysis_sets <- function(surv, outcomes, params) {

  sets <- vector("list", nrow(outcomes))
  names(sets) <- outcomes$key

  for (i in seq_len(nrow(outcomes))) {
    o <- outcomes[i, ]
    prev_cols <- .resolve_prev(o, outcomes)
    d <- exclude_prevalent_clusters(surv, prev_cols, params)

    cr <- d[[paste0(o$key, "_incid_cr")]]
    tt <- d[[paste0("time_to_", o$key)]]
    sizes <- table(d[["subclass"]])

    sets[[i]] <- structure(
      list(data = d,
           key = o$key,
           prev_cols = prev_cols,
           n = nrow(d),
           n_clusters = length(unique(d[["subclass"]])),
           cases = sum(!is.na(cr) & cr == 1),
           competing = sum(!is.na(cr) & cr == 2),
           n_missing_time = sum(is.na(tt)),
           cluster_size = as.integer(sizes)),
      class = "tte_analysis_set")
  }

  sets
}
#' @noRd
.resolve_prev <- function(o, outcomes) {
  if (!is.na(o$prev[[1L]])) return(o$prev[[1L]])
  comps <- o$components[[1L]]
  if (is.null(comps)) return(character(0))
  unique(stats::na.omit(outcomes$prev[match(comps, outcomes$key)]))
}
#' @noRd
exclude_prevalent_clusters <- function(data, prev_cols, params) {

  if (!length(prev_cols)) return(data)
  prev_cols <- intersect(prev_cols, names(data))
  if (!length(prev_cols)) return(data)

  affected <- rep(FALSE, nrow(data))
  for (cl in prev_cols) {
    v <- data[[cl]]
    flag <- if (is.factor(v)) {
      !is.na(v) & as.character(v) == "1"
    } else {
      !is.na(v) & v == 1
    }
    affected <- affected | flag
  }

  bad_clusters <- unique(data[["subclass"]][affected])
  drop_idx <- which(data[["subclass"]] %in% bad_clusters)
  drop_rows(data, drop_idx)
}
#' @noRd
tte_fit_outcomes <- function(sets, outcomes, params, strict = TRUE) {

  fits <- vector("list", nrow(outcomes))
  names(fits) <- outcomes$key

  for (i in seq_len(nrow(outcomes))) {
    o <- outcomes[i, ]
    fits[[i]] <- fit_outcome(sets[[o$key]], o, params, strict = strict)
  }

  fits
}
#' @noRd
fit_outcome <- function(set, outcome, params, strict = TRUE) {

  key <- outcome$key[[1L]]
  model <- outcome$model[[1L]]
  fail <- function(reason) {
    structure(list(ok = FALSE, reason = reason, key = key), class = "tte_fit")
  }

  if (is.null(set) || !set$n) return(fail("empty analysis set"))
  if (set$cases < params$min_events) {
    return(fail(sprintf("insufficient events: %d < min_events = %g",
                        set$cases, params$min_events)))
  }

  d <- set$data
  tcol <- paste0("time_to_", key)
  ccol <- paste0(key, "_incid_cr")
  if (!all(c(tcol, ccol) %in% names(d))) {
    return(fail(paste0("derived columns missing for \"", key, "\"")))
  }

  keep <- !is.na(d[[tcol]]) & !is.na(d[[ccol]]) & d[[tcol]] > 0
  n_dropped <- sum(!keep)
  if (n_dropped > 0 && n_dropped / length(keep) > 0.05) {
    warning("Outcome \"", key, "\": dropping ", n_dropped, " of ",
            length(keep), " rows with missing or non-positive follow-up ",
            "time. If `followup_end` is NULL, censored participants whose ",
            "date column is NA get NA follow-up and are silently excluded ",
            "-- for mortality this removes everyone still alive. Set ",
            "`followup_end` to the administrative censoring date or column.",
            call. = FALSE)
  }
  d <- d[keep, ]
  if (!nrow(d)) return(fail("no rows with positive, non-missing follow-up"))

  out <- tryCatch({
    if (identical(model, "cox")) {
      .fit_cox(d, tcol, ccol, params)
    } else {
      .fit_crr(d, tcol, ccol, params)
    }
  }, error = function(e) {
    if (isTRUE(strict)) {
      stop("Outcome \"", key, "\" failed: ", conditionMessage(e),
           "\n(Set strict = FALSE to record the failure and continue.)",
           call. = FALSE)
    }
    fail(conditionMessage(e))
  })

  out$key <- key
  out
}
#' @noRd
.fit_crr <- function(d, tcol, ccol, params) {
  cov1 <- as.matrix(data.frame(exposure = as.numeric(d[[params$exposure_col]])))
  fit <- crrSC::crrc(ftime = as.numeric(d[[tcol]]),
                     fstatus = as.integer(d[[ccol]]),
                     cov1 = cov1,
                     cluster = as.integer(factor(d[["subclass"]])))
  b <- as.numeric(fit$coef)[1L]
  se <- sqrt(as.numeric(diag(as.matrix(fit$var)))[1L])
  z <- b / se
  structure(list(ok = TRUE,
                 est = exp(b),
                 lower = exp(b - stats::qnorm(0.975) * se),
                 upper = exp(b + stats::qnorm(0.975) * se),
                 p = 2 * stats::pnorm(-abs(z)),
                 coef = b, se = se, model = "crr", fit = fit),
            class = "tte_fit")
}
#' @noRd
.fit_cox <- function(d, tcol, ccol, params) {
  ev <- as.integer(d[[ccol]] == 1)
  df <- data.frame(time = as.numeric(d[[tcol]]),
                   event = ev,
                   exposure = as.numeric(d[[params$exposure_col]]),
                   subclass = d[["subclass"]])
  ## The cluster= argument, not a cluster() term in the formula: coxph detects
  ## cluster() as a formula *special* by name, and the namespace-qualified form
  ## survival::cluster() is not matched by that detection, so clustered robust
  ## SEs would silently not be applied as intended.
  fit <- survival::coxph(
    survival::Surv(time, event) ~ exposure,
    data = df, cluster = subclass)
  s <- summary(fit)
  b <- unname(s$coefficients["exposure", "coef"])
  se_col <- if ("robust se" %in% colnames(s$coefficients)) "robust se" else "se(coef)"
  se <- unname(s$coefficients["exposure", se_col])
  z <- b / se
  structure(list(ok = TRUE,
                 est = exp(b),
                 lower = exp(b - stats::qnorm(0.975) * se),
                 upper = exp(b + stats::qnorm(0.975) * se),
                 p = 2 * stats::pnorm(-abs(z)),
                 coef = b, se = se, model = "cox", fit = fit),
            class = "tte_fit")
}
#' @noRd
tte_assemble_results <- function(fits, sets, outcomes, params,
                                 sort_by = c("est", "registry")) {

  sort_by <- match.arg(sort_by)

  rows <- lapply(seq_len(nrow(outcomes)), function(i) {
    o <- outcomes[i, ]
    summarise_outcome(fits[[o$key]], sets[[o$key]], o)
  })
  res <- do.call(rbind, rows)

  res$q <- NA_real_
  ok <- !is.na(res$p)
  if (any(ok)) res$q[ok] <- stats::p.adjust(res$p[ok], method = "fdr")

  if (identical(sort_by, "est")) {
    res <- res[order(res$est, na.last = TRUE), , drop = FALSE]
  }
  rownames(res) <- NULL

  class(res) <- unique(c("tte_results", class(res)))
  res
}
#' Run the pipeline
#'
#' Runs every arm through the six stages and returns all intermediates, so a
#' question like "did matching go wrong or did the model?" is answerable
#' without re-running anything.
#'
#' Arms whose source file and matching-relevant parameters hash to the same
#' cache key share one matched set, computed once. That is what makes an
#' intention-to-treat arm and a per-protocol arm built from the same file
#' *guaranteed* to have identical pairs.
#'
#' @param scenarios A [tte_scenarios()] object.
#' @param outcomes An outcome registry from [tte_outcomes()].
#' @param params A [tte_params()] object.
#' @param cache_dir Directory for matched-set checkpoints, or `NULL`. Keys
#'   cover the source file\'s content signature and the matching parameters, so
#'   an edited file or a changed caliper cannot silently return a stale set.
#' @param out_dir Directory for `<label>_surv.rds` and `<label>_results.csv`,
#'   or `NULL` to write nothing.
#' @param through Stop after a given stage and return what exists so far:
#'   `"match"`, `"surv"`, `"sets"`, `"fits"`, or `"results"` (default). Use it
#'   to inspect an intermediate without paying for the rest.
#' @param matched Optionally inject an already-constructed matched set,
#'   bypassing stages 1 and 2 — the route for intervening between stages.
#'
#' @return A list with `arms` (per-arm `matched`, `surv`, `sets`, `fits`,
#'   `results`), `results` (the tables alone), `reference` and `through`.
#' @examples
#' set.seed(1); n <- 200
#' cohort <- data.frame(
#'   pid = sprintf("P%03d", seq_len(n)),
#'   exposed = rep(c(1L, 0L), c(60, 140)),
#'   age = round(stats::rnorm(n, 60, 8)), biom = stats::rnorm(n, 6, 0.6),
#'   base_date = as.Date("2010-01-01"),
#'   dth = stats::rbinom(n, 1, 0.2), dth_date = as.Date("2015-06-01"),
#'   oa_prev = stats::rbinom(n, 1, 0.1), oa_status = 0L,
#'   oa_incid = stats::rbinom(n, 1, 0.3), oa_date = as.Date("2014-03-01"))
#' outcomes <- tte_outcomes("oa", "Osteoarthritis", "oa_prev", "oa_status",
#'                          "oa_incid", "oa_date", model = "cox")
#' params <- tte_params("pid", "exposed", exposed ~ age + biom, "base_date",
#'                      death_status_col = "dth", death_date_col = "dth_date",
#'                      match_calipers = c(biom = 0.5), ps_caliper = 0.2,
#'                      followup_end = as.Date("2019-12-31"))
#'
#' f <- tempfile(fileext = ".csv")
#' utils::write.csv(cohort, f, row.names = FALSE)
#' scen <- tte_scenarios(c(A = f, B = f),
#'                       list(tte_strategy("itt"), tte_strategy("itt", label = "ITT2")))
#'
#' run <- tte_run(scen, outcomes, params)
#' run$results$A
#'
#' ## both arms read one file, so the pairs are identical by construction
#' identical(run$arms$A$matched$user_ids, run$arms$B$matched$user_ids)
#'
#' ## stop early to inspect the matched set alone
#' m <- tte_run(scen, outcomes, params, through = "match")
#' m$arms$A$matched
#' unlink(f)
#' @export
tte_run <- function(scenarios, outcomes, params, cache_dir = NULL,
                    out_dir = NULL,
                    through = c("results", "match", "surv", "sets", "fits"),
                    matched = NULL) {

  through <- match.arg(through)

  if (!inherits(scenarios, "tte_scenarios")) {
    stop("`scenarios` must come from tte_scenarios().", call. = FALSE)
  }

  keys <- vapply(scenarios$source, function(s) {
    if (is.character(s) && length(s) == 1L) tte_cache_key(s, params) else NA_character_
  }, character(1))
  shared <- list()
  arms <- list()

  for (i in seq_along(scenarios$label)) {
    k <- keys[i]
    m <- if (!is.na(k)) shared[[k]] else NULL
    arms[[scenarios$label[i]]] <- .run_scenario(
      label = scenarios$label[i],
      source = scenarios$source[i],
      strategy = scenarios$strategy[[i]],
      outcomes = outcomes,
      params = params,
      cache_dir = cache_dir,
      out_dir = out_dir,
      matched = if (!is.null(matched)) matched else m,
      through = through)
    if (is.null(m) && !is.na(k)) shared[[k]] <- arms[[scenarios$label[i]]]$matched
  }

  list(arms = arms,
       results = lapply(arms, `[[`, "results"),
       reference = scenarios$reference,
       through = through)
}
#' @noRd
.run_scenario <- function(label, source, strategy, outcomes, params,
                          cache_dir = NULL, out_dir = NULL,
                          matched = NULL, through = "results") {

  message("[", label, "] stage 1: cohort")
  cohort <- .load_cohort(source, outcomes, params)

  if (is.null(matched)) {
    message("[", label, "] stage 2: matching")
    matched <- tte_match(cohort, params, cache_dir = cache_dir,
                         source = source)
  } else {
    message("[", label, "] stage 2: reusing shared matched set")
  }

  if (identical(through, "match")) {
    return(list(label = label, matched = matched, surv = NULL, sets = NULL,
                fits = NULL, results = NULL))
  }

  message("[", label, "] stage 3: analysis time (", strategy$label, ")")
  surv <- tte_analysis_time(matched, outcomes, strategy, params)
  attr(surv, "params") <- params

  if (identical(through, "surv")) {
    return(list(label = label, matched = matched, surv = surv, sets = NULL,
                fits = NULL, results = NULL))
  }

  message("[", label, "] stage 4: per-outcome sets")
  sets <- tte_analysis_sets(surv, outcomes, params)

  if (identical(through, "sets")) {
    return(list(label = label, matched = matched, surv = surv, sets = sets,
                fits = NULL, results = NULL))
  }

  message("[", label, "] stage 5: fits")
  fits <- tte_fit_outcomes(sets, outcomes, params)

  if (identical(through, "fits")) {
    return(list(label = label, matched = matched, surv = surv, sets = sets,
                fits = fits, results = NULL))
  }

  message("[", label, "] stage 6: results")
  results <- tte_assemble_results(fits, sets, outcomes, params)

  if (!is.null(out_dir)) {
    if (!dir.exists(out_dir)) {
      dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    }
    slug <- gsub("[^A-Za-z0-9]+", "_", label)
    saveRDS(surv, file.path(out_dir, paste0(slug, "_surv.rds")))
    utils::write.csv(as.data.frame(results),
                     file.path(out_dir, paste0(slug, "_results.csv")),
                     row.names = FALSE)
  }

  list(label = label, matched = matched, surv = surv, sets = sets,
       fits = fits, results = results)
}
