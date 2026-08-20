# ---------------------------------------------------------------------------
# Output and inspection: plots, diagnostics, the escape hatch, and the
# baseline comparison used during a migration.
# ---------------------------------------------------------------------------

#' Forest plot
#'
#' Draws one arm, or several side by side when given a named list of results
#' tables.
#'
#' Arms are merged on the outcome **key**, never on row position: a divergence
#' in row order between arms would otherwise relabel every row rather than
#' raise an error. Row order is recomputed from the reference arm on each run,
#' so it cannot go stale. Outcomes with `in_plot = FALSE` are excluded.
#'
#' @param results A results table from [tte_run()], or a named list of them.
#' @param scenarios Optional [tte_scenarios()] object supplying the reference
#'   arm.
#' @param reference Arm label whose ordering governs. Defaults to
#'   `scenarios$reference`, else the first arm.
#' @param title Plot title.
#' @param ref_line Reference line. Default 1.
#'
#' @return A `forestploter` object.
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
#' f <- tempfile(fileext = ".csv")
#' utils::write.csv(cohort, f, row.names = FALSE)
#' scen <- tte_scenarios(c(A = f), list(tte_strategy("itt")))
#' run <- tte_run(scen, outcomes, params)
#'
#' p <- tte_forest(run$results$A, title = "Demo")
#' class(p)
#'
#' ## several arms: merged by outcome key, ordered by the reference arm
#' q <- tte_forest(run$results, scen)
#' class(q)
#' unlink(f)
#' @export
tte_forest <- function(results, scenarios = NULL, reference = NULL,
                       title = NULL, ref_line = 1) {
  if (inherits(results, "tte_results") || is.data.frame(results)) {
    return(.forest_single(results, title, ref_line))
  }
  if (!is.list(results) || !length(results)) {
    stop("`results` must be a results table or a named list of them.",
         call. = FALSE)
  }
  .forest_combined(results, scenarios, reference, title, ref_line)
}
#' @noRd
.forest_single <- function(results, title, ref_line) {

  d <- results[isTRUE_vec(results$in_plot), , drop = FALSE]
  if (!nrow(d)) stop("No outcomes have in_plot = TRUE.", call. = FALSE)

  tab <- data.frame(
    Outcome = d$label,
    N = d$n,
    Cases = d$cases,
    ` ` = paste(rep(" ", 24), collapse = ""),
    `HR (95% CI)` = ifelse(is.na(d$ci), "not estimated", d$ci),
    p = .fmt_p(d$p),
    q = .fmt_p(d$q),
    check.names = FALSE, stringsAsFactors = FALSE)

  forestploter::forest(tab,
                       est = d$est, lower = d$lower, upper = d$upper,
                       ci_column = which(names(tab) == " "),
                       ref_line = ref_line,
                       title = title)
}
#' @noRd
.forest_combined <- function(results_list, scenarios, reference, title,
                             ref_line) {

  if (is.null(reference)) {
    reference <- if (!is.null(scenarios)) scenarios$reference
                 else names(results_list)[1L]
  }
  if (!reference %in% names(results_list)) {
    stop("Reference arm \"", reference, "\" is not among: ",
         paste(names(results_list), collapse = ", "), call. = FALSE)
  }

  ref <- results_list[[reference]]
  ref <- ref[isTRUE_vec(ref$in_plot), , drop = FALSE]
  key_order <- ref$key[order(ref$est, na.last = TRUE)]

  tab <- data.frame(Outcome = ref$label[match(key_order, ref$key)],
                    stringsAsFactors = FALSE, check.names = FALSE)
  est <- lower <- upper <- list()
  ci_cols <- integer(0)

  for (arm in names(results_list)) {
    r <- results_list[[arm]]
    idx <- match(key_order, r$key)
    if (anyNA(idx)) {
      warning("Arm \"", arm, "\" is missing outcome(s): ",
              paste(key_order[is.na(idx)], collapse = ", "), call. = FALSE)
    }
    tab[[paste0(arm, " ")]] <- paste(rep(" ", 24), collapse = "")
    ci_cols <- c(ci_cols, ncol(tab))
    tab[[paste0(arm, " HR (95% CI)")]] <-
      ifelse(is.na(r$ci[idx]), "not estimated", r$ci[idx])
    est[[arm]] <- r$est[idx]
    lower[[arm]] <- r$lower[idx]
    upper[[arm]] <- r$upper[idx]
  }

  forestploter::forest(tab,
                       est = est, lower = lower, upper = upper,
                       ci_column = ci_cols,
                       ref_line = ref_line,
                       title = title)
}
#' Balance diagnostics
#'
#' A non-blocking branch off matching: nothing downstream depends on it.
#' Reports the match tally, standardised mean differences before and after
#' matching, and a pairing scatter for each caliper covariate — generalised to
#' whatever `match_calipers` names, so no covariate is privileged in code.
#'
#' Standardised differences are computed internally: pooled-SD for continuous
#' variables, the Yang & Dalton multi-category formulation for factors — the
#' same definitions `tableone` uses, so values stay comparable with SMD tables
#' from an original pipeline.
#'
#' @param matched A matched set, e.g. `run$arms$A$matched`.
#' @param cohort The pre-matching cohort, for the "before" comparison.
#' @param params A [tte_params()] object.
#' @param out_dir Optional directory; scatter plots are written there as PNG.
#'
#' @return Invisibly, a list with `tally`, `smd_after`, `smd_before` and
#'   `scatter_data`.
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
#' f <- tempfile(fileext = ".csv")
#' utils::write.csv(cohort, f, row.names = FALSE)
#' scen <- tte_scenarios(c(A = f), list(tte_strategy("itt")))
#' run <- tte_run(scen, outcomes, params)
#'
#' bal <- tte_balance_report(run$arms$A$matched, params = params)
#' bal$tally
#' bal$smd_after
#' unlink(f)
#' @export
tte_balance_report <- function(matched, cohort = NULL, params,
                               out_dir = NULL) {

  d <- matched$data
  exp_col <- params$exposure_col

  tally <- table(matched$match_num)

  vars <- setdiff(all.vars(params$ps_formula), exp_col)
  vars <- unique(c(vars, names(params$match_calipers)))

  smd_after <- .smd_table(d, vars, exp_col)
  smd_before <- if (!is.null(cohort)) {
    .smd_table(cohort, intersect(vars, names(cohort)), exp_col)
  } else NULL

  scatter <- list()
  for (cl in names(params$match_calipers)) {
    if (!cl %in% names(d)) next
    ex <- !is.na(d[[exp_col]]) & d[[exp_col]] == 1
    u <- stats::setNames(as.numeric(d[[cl]])[ex],
                         as.character(d[["subclass"]][ex]))
    c_ <- stats::setNames(as.numeric(d[[cl]])[!ex],
                          as.character(d[["subclass"]][!ex]))
    common <- intersect(names(u), names(c_))
    scatter[[cl]] <- data.frame(subclass = common,
                                exposed = unname(u[common]),
                                control = unname(c_[common]),
                                stringsAsFactors = FALSE)
    if (!is.null(out_dir)) {
      if (!dir.exists(out_dir)) {
        dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
      }
      grDevices::png(file.path(out_dir, paste0("pairing_", cl, ".png")),
                     width = 800, height = 800)
      on.exit(grDevices::dev.off(), add = TRUE)
      graphics::plot(scatter[[cl]]$exposed, scatter[[cl]]$control,
                     xlab = paste(cl, "(exposed)"),
                     ylab = paste(cl, "(control)"),
                     main = paste("Pairing on", cl), pch = 16,
                     col = grDevices::rgb(0, 0, 0, 0.3))
      graphics::abline(0, 1, col = "red")
    }
  }

  invisible(list(tally = tally, smd_after = smd_after,
                 smd_before = smd_before, scatter_data = scatter))
}
#' Processed analysis data
#'
#' The escape hatch. Returns the derived dataset so analyses this package does
#' not perform — inverse-probability weighting or censoring weights in
#' particular — can be carried out downstream.
#'
#' The attributes carry what the data frame alone cannot tell you: which
#' estimand the event times encode, that pairs are **not** independent
#' observations and any model must account for `subclass`, the time-zero rule
#' in words, how administrative censoring was applied, and for per-protocol
#' arms the per-participant dispensing summary — the natural covariate set for
#' modelling discontinuation.
#'
#' @param x An arm from [tte_run()], e.g. `run$arms$A`.
#' @param outcome `NULL` (default) for the full derived dataset, or an outcome
#'   key for that outcome\'s post-exclusion set. The populations differ: pairs
#'   affected at baseline are dropped per outcome.
#'
#' @return A `data.frame` with attributes as described.
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
#' f <- tempfile(fileext = ".csv")
#' utils::write.csv(cohort, f, row.names = FALSE)
#' scen <- tte_scenarios(c(A = f), list(tte_strategy("itt")))
#' run <- tte_run(scen, outcomes, params)
#'
#' d <- tte_analysis_data(run$arms$A)
#' attr(d, "cluster_col")
#' attr(d, "time_zero")
#' unlink(f)
#' @export
tte_analysis_data <- function(x, outcome = NULL) {

  surv <- if (inherits(x, "tte_surv")) x else x$surv
  sets <- if (inherits(x, "tte_surv")) NULL else x$sets
  params <- attr(surv, "params")

  d <- if (is.null(outcome)) {
    surv
  } else {
    if (is.null(sets) || is.null(sets[[outcome]])) {
      stop("No analysis set for outcome \"", outcome,
           "\". Pass an arm from tte_run().",
           call. = FALSE)
    }
    sets[[outcome]]$data
  }

  d <- as.data.frame(d)
  attr(d, "strategy") <- attr(surv, "strategy_label")
  attr(d, "id_col") <- if (!is.null(params)) params$id_col else NA_character_
  attr(d, "cluster_col") <- "subclass"
  attr(d, "time_zero") <- attr(surv, "time_zero_rule")
  attr(d, "followup_end") <- attr(surv, "followup_end")
  attr(d, "rx_summary") <- attr(surv, "rx_summary")
  attr(d, "outcome") <- outcome
  d
}
#' Compare a run against a saved baseline
#'
#' Comparison is layered so a difference can be attributed to a stage.
#' False-discovery-rate-adjusted p-values are their own layer: correcting a
#' single outcome shifts every q-value as an arithmetic consequence, and that
#' must read as expected propagation rather than as a second, unexplained
#' difference.
#'
#' Layers absent from either side are skipped.
#'
#' @param current,baseline Lists with any of `pairs` (a two-column
#'   cluster/id mapping), `ps` (named numeric), `counts` (data.frame keyed by
#'   outcome), `estimates` (with `key`, `est`, `lower`, `upper`, `p`), `q`
#'   (with `key`, `q`), `row_order` (character vector of keys).
#' @param tol Tolerance for the floating-point layers.
#'
#' @return An object of class `tte_comparison`.
#' @examples
#' base <- list(
#'   estimates = data.frame(key = c("a", "b"), est = c(1.10, 0.90),
#'                          lower = c(0.9, 0.8), upper = c(1.3, 1.0),
#'                          p = c(0.20, 0.04)),
#'   q = data.frame(key = c("a", "b"), q = c(0.20, 0.08)),
#'   row_order = c("b", "a"))
#'
#' tte_compare(base, base)
#'
#' ## one corrected estimate shifts every q-value; the layers keep that
#' ## expected propagation separate from the real change
#' cur <- base
#' cur$estimates$est[1] <- 1.35
#' cur$q$q <- c(0.15, 0.06)
#' tte_compare(cur, base)
#' @export
tte_compare <- function(current, baseline, tol = 1e-8) {

  layer <- function(name, pass, detail) {
    list(name = name, pass = pass, detail = detail)
  }
  out <- list()

  if (!is.null(current$pairs) && !is.null(baseline$pairs)) {
    a <- do.call(paste, c(as.list(as.data.frame(current$pairs)), sep = "|"))
    b <- do.call(paste, c(as.list(as.data.frame(baseline$pairs)), sep = "|"))
    only_cur <- setdiff(a, b); only_base <- setdiff(b, a)
    out$pairs <- layer(
      "matched pair IDs (exact)",
      length(only_cur) == 0L && length(only_base) == 0L,
      sprintf("%d only in current, %d only in baseline",
              length(only_cur), length(only_base)))
  }

  num_layer <- function(cur, base, nm) {
    common <- intersect(names(cur), names(base))
    if (!length(common)) return(layer(nm, FALSE, "no common identifiers"))
    d <- abs(cur[common] - base[common])
    layer(nm, all(d <= tol, na.rm = TRUE) &&
            length(common) == length(base),
          sprintf("max abs diff %.3g over %d values (%d missing)",
                  suppressWarnings(max(d, na.rm = TRUE)), length(common),
                  length(base) - length(common)))
  }

  if (!is.null(current$ps) && !is.null(baseline$ps)) {
    out$ps <- num_layer(current$ps, baseline$ps,
                        sprintf("propensity scores (tol %.0e)", tol))
  }

  if (!is.null(current$counts) && !is.null(baseline$counts)) {
    cc <- current$counts[order(current$counts$key), , drop = FALSE]
    bb <- baseline$counts[order(baseline$counts$key), , drop = FALSE]
    same <- isTRUE(all.equal(cc, bb, check.attributes = FALSE))
    out$counts <- layer("per-outcome counts (exact)", same,
                        if (same) "identical" else "see all.equal() output")
  }

  if (!is.null(current$estimates) && !is.null(baseline$estimates)) {
    cols <- intersect(c("est", "lower", "upper", "p"),
                      intersect(names(current$estimates),
                                names(baseline$estimates)))
    worst <- 0; ok <- TRUE
    for (cl in cols) {
      cur <- stats::setNames(current$estimates[[cl]], current$estimates$key)
      base <- stats::setNames(baseline$estimates[[cl]], baseline$estimates$key)
      common <- intersect(names(cur), names(base))
      d <- suppressWarnings(max(abs(cur[common] - base[common]), na.rm = TRUE))
      if (is.finite(d)) worst <- max(worst, d)
      ok <- ok && isTRUE(d <= tol)
    }
    out$estimates <- layer(
      sprintf("estimates and raw p (tol %.0e)", tol), ok,
      sprintf("max abs diff %.3g across %s", worst, paste(cols, collapse = "/")))
  }

  if (!is.null(current$q) && !is.null(baseline$q)) {
    cur <- stats::setNames(current$q$q, current$q$key)
    base <- stats::setNames(baseline$q$q, baseline$q$key)
    out$q <- num_layer(cur, base, sprintf("FDR q-values (tol %.0e)", tol))
  }

  if (!is.null(current$row_order) && !is.null(baseline$row_order)) {
    same <- identical(as.character(current$row_order),
                      as.character(baseline$row_order))
    out$row_order <- layer("forest row order (exact)", same,
                           if (same) "identical" else
                             paste("current:",
                                   paste(current$row_order, collapse = ",")))
  }

  structure(list(layers = out, tol = tol), class = "tte_comparison")
}
#' @noRd
drop_rows <- function(x, i) {
  i <- i[!is.na(i)]
  if (!length(i)) return(x)
  if (is.data.frame(x)) x[-i, , drop = FALSE] else x[-i]
}
#' @noRd
.as_date <- function(x) {
  if (inherits(x, "Date")) return(x)
  if (inherits(x, "POSIXt")) return(as.Date(x))
  suppressWarnings(as.Date(as.character(x)))
}
#' @noRd
require_cols <- function(data, outcomes, params, error = TRUE) {
  need <- c(params$id_col, params$exposure_col,
            params$time_zero_col,
            params$death_status_col, params$death_date_col,
            all.vars(params$ps_formula),
            names(params$match_calipers))
  if (is.character(params$followup_end)) need <- c(need, params$followup_end)

  elem <- outcomes[outcomes$model != "derived", ]
  need <- c(need, elem$prev, elem$status, elem$incid, elem$date)

  der <- outcomes[outcomes$model == "derived", ]
  if (nrow(der)) {
    need <- c(need, unlist(der$dates), stats::na.omit(der$prev))
  }

  need <- unique(stats::na.omit(need))
  missing_cols <- setdiff(need, names(data))
  if (length(missing_cols)) {
    msg <- paste0("Missing required column(s): ",
                  paste(missing_cols, collapse = ", "))
    if (error) stop(msg, call. = FALSE) else warning(msg, call. = FALSE)
  }
  invisible(missing_cols)
}
#' @noRd
.matching_params <- function(params) {
  list(ps_formula = paste(deparse(params$ps_formula), collapse = " "),
       match_calipers = params$match_calipers,
       ps_caliper = params$ps_caliper,
       id_col = params$id_col,
       exposure_col = params$exposure_col)
}
#' @noRd
tte_cache_key <- function(source, params) {
  src <- normalizePath(source, mustWork = FALSE)
  info <- file.info(src)
  digest::digest(list(source = src,
                      size = info$size,
                      mtime = info$mtime,
                      params = .matching_params(params)))
}
#' @noRd
tte_cache <- function(source, params, cache_dir, compute, force = FALSE) {
  if (is.null(cache_dir)) return(compute())

  key <- tte_cache_key(source, params)
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  }
  path <- file.path(cache_dir, paste0("matched_", key, ".rds"))

  if (!force && file.exists(path)) {
    obj <- readRDS(path)
    if (!identical(obj$cache_key, key)) {
      stop("Cache file ", path, " has key ", obj$cache_key,
           " but the current configuration hashes to ", key,
           ". Delete the file or call with force = TRUE.", call. = FALSE)
    }
    if (!isTRUE(all.equal(obj$params_used, .matching_params(params)))) {
      stop("Cached matched set was built with different matching parameters. ",
           "Delete ", path, " or call with force = TRUE.", call. = FALSE)
    }
    return(obj)
  }

  obj <- compute()
  obj$cache_key <- key
  obj$params_used <- .matching_params(params)
  saveRDS(obj, path)
  obj
}
#' @noRd
derive_outcome <- function(data, key, incid, date, year_days,
                           followup_end = NULL,
                           death_status = NULL, death_date = NULL,
                           censor_date = NULL) {

  t0 <- .as_date(data[["time_zero"]])
  ev <- !is.na(data[[incid]]) & data[[incid]] == 1
  ev_date <- .as_date(data[[date]])
  ev[is.na(ev_date)] <- FALSE

  dead <- rep(FALSE, nrow(data))
  d_date <- rep(as.Date(NA), nrow(data))
  if (!is.null(death_status) && !is.null(death_date)) {
    dead <- !is.na(data[[death_status]]) & data[[death_status]] == 1
    d_date <- .as_date(data[[death_date]])
    dead[is.na(d_date)] <- FALSE
  }

  end_date <- if (is.null(followup_end)) {
    ev_date
  } else if (inherits(followup_end, "Date")) {
    rep(followup_end, nrow(data))
  } else {
    .as_date(data[[followup_end]])
  }

  cr <- integer(nrow(data))
  cr[!ev & dead] <- 2L
  cr[ev] <- 1L

  end <- end_date
  end[!ev & dead] <- d_date[!ev & dead]
  end[ev] <- ev_date[ev]

  ## Cap at the strategy's censoring date. Events and competing deaths after
  ## that date have already been removed by apply_censoring(); without this cap
  ## those participants would be carried on to administrative end and gain
  ## person-time they never contributed.
  if (!is.null(censor_date)) {
    cd <- .as_date(censor_date)
    late <- !is.na(cd) & !is.na(end) & end > cd
    end[late] <- cd[late]
  }

  data[[paste0("time_to_", key)]] <-
    as.numeric(end - t0) / year_days
  data[[paste0(key, "_incid_cr")]] <- cr
  data
}
#' @noRd
apply_censoring <- function(data, key, status, incid, date, censor_date) {
  d <- .as_date(data[[date]])
  cd <- .as_date(censor_date)
  after <- !is.na(d) & !is.na(cd) & d > cd
  if (any(after)) {
    if (!is.na(incid) && incid %in% names(data)) data[[incid]][after] <- 0
    if (!is.na(status) && status %in% names(data)) data[[status]][after] <- 0
    if (!is.na(date) && date %in% names(data)) data[[date]][after] <- NA
  }
  data
}
#' @noRd
derive_composite_endpoint <- function(data, spec, registry, params) {

  key <- spec$key[[1L]]
  comps <- spec$components[[1L]]
  from <- spec$from[[1L]]
  dates <- spec$dates[[1L]]

  miss <- setdiff(comps, registry$key)
  if (length(miss)) {
    stop("Derived outcome \"", key, "\" names unknown component(s): ",
         paste(miss, collapse = ", "), call. = FALSE)
  }
  if (is.null(dates)) {
    dates <- registry$date[match(comps, registry$key)]
  }

  ev <- rep(FALSE, nrow(data))
  for (i in seq_along(comps)) {
    col <- if (identical(from[i], "cr")) {
      paste0(comps[i], "_incid_cr")
    } else {
      registry$incid[match(comps[i], registry$key)]
    }
    if (!col %in% names(data)) {
      stop("Derived outcome \"", key, "\" needs column \"", col,
           "\", which is not present.", call. = FALSE)
    }
    ev <- ev | (!is.na(data[[col]]) & data[[col]] == 1)
  }

  dmat <- do.call(cbind, lapply(dates, function(cl) as.numeric(.as_date(data[[cl]]))))
  first <- suppressWarnings(apply(dmat, 1L, function(r) {
    r <- r[!is.na(r)]
    if (!length(r)) NA_real_ else min(r)
  }))
  first <- as.Date(first, origin = "1970-01-01")

  t0 <- .as_date(data[["time_zero"]])
  dead <- rep(FALSE, nrow(data)); d_date <- rep(as.Date(NA), nrow(data))
  if (!is.null(params$death_status_col) && !is.null(params$death_date_col)) {
    dead <- !is.na(data[[params$death_status_col]]) &
      data[[params$death_status_col]] == 1
    d_date <- .as_date(data[[params$death_date_col]])
    dead[is.na(d_date)] <- FALSE
  }

  ev <- ev & !is.na(first)
  cr <- integer(nrow(data))
  cr[!ev & dead] <- 2L
  cr[ev] <- 1L

  end <- first
  end[!ev & dead] <- d_date[!ev & dead]
  if (!is.null(params$followup_end)) {
    fe <- if (inherits(params$followup_end, "Date")) {
      rep(params$followup_end, nrow(data))
    } else .as_date(data[[params$followup_end]])
    idx <- !ev & !dead
    end[idx] <- fe[idx]
  }

  ## same censoring cap as derive_outcome(); see its documentation
  if ("censor_date" %in% names(data)) {
    cd <- .as_date(data[["censor_date"]])
    late <- !is.na(cd) & !is.na(end) & end > cd
    end[late] <- cd[late]
  }

  data[[paste0("time_to_", key)]] <- as.numeric(end - t0) / params$year_days
  data[[paste0(key, "_incid_cr")]] <- cr
  data
}
#' @noRd
summarise_outcome <- function(fit, set, outcome, digits = 2) {
  data.frame(
    key = outcome$key[[1L]],
    label = outcome$label[[1L]],
    n = set$n,
    n_clusters = set$n_clusters,
    cases = set$cases,
    competing = set$competing,
    est = if (isTRUE(fit$ok)) fit$est else NA_real_,
    lower = if (isTRUE(fit$ok)) fit$lower else NA_real_,
    upper = if (isTRUE(fit$ok)) fit$upper else NA_real_,
    p = if (isTRUE(fit$ok)) fit$p else NA_real_,
    ci = if (isTRUE(fit$ok)) {
      format_ci(fit$est, fit$lower, fit$upper, digits)
    } else NA_character_,
    ok = isTRUE(fit$ok),
    reason = if (isTRUE(fit$ok)) NA_character_ else fit$reason,
    in_plot = outcome$in_plot[[1L]],
    stringsAsFactors = FALSE
  )
}
#' @noRd
format_ci <- function(est, lower, upper, digits = 2) {
  f <- function(x) formatC(x, format = "f", digits = digits)
  paste0(f(est), " (", f(lower), "-", f(upper), ")")
}
#' @noRd
.smd_table <- function(data, vars, exp_col) {
  vars <- intersect(vars, names(data))
  vars <- setdiff(vars, exp_col)
  if (!length(vars)) return(NULL)

  data <- as.data.frame(data)
  g <- !is.na(data[[exp_col]]) & data[[exp_col]] == 1

  smd <- vapply(vars, function(v) {
    x <- data[[v]]
    if (is.numeric(x)) .smd_continuous(x[g], x[!g]) else
      .smd_categorical(x[g], x[!g])
  }, numeric(1))

  data.frame(variable = vars, smd = unname(smd),
             row.names = NULL, stringsAsFactors = FALSE)
}
#' @noRd
.smd_continuous <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  if (!length(a) || !length(b)) return(NA_real_)
  s <- sqrt((stats::var(a) + stats::var(b)) / 2)
  if (!is.finite(s) || s == 0) return(0)
  abs(mean(a) - mean(b)) / s
}
#' @noRd
.smd_categorical <- function(a, b) {
  a <- as.character(a[!is.na(a)]); b <- as.character(b[!is.na(b)])
  if (!length(a) || !length(b)) return(NA_real_)
  lev <- sort(unique(c(a, b)))
  if (length(lev) < 2L) return(0)
  lev <- lev[-1L]                      # drop one level: K-1 dimensions
  p1 <- vapply(lev, function(l) mean(a == l), numeric(1))
  p2 <- vapply(lev, function(l) mean(b == l), numeric(1))
  d <- p1 - p2
  cov_of <- function(p) {
    s <- -outer(p, p)
    diag(s) <- p * (1 - p)
    s
  }
  S <- (cov_of(p1) + cov_of(p2)) / 2
  out <- tryCatch(sqrt(drop(t(d) %*% solve(S) %*% d)),
                  error = function(e) NA_real_)
  if (!is.finite(out)) NA_real_ else out
}
#' @noRd
isTRUE_vec <- function(x) !is.na(x) & x

#' @keywords internal
.fmt_p <- function(p) {
  ifelse(is.na(p), "",
         ifelse(p < 0.001, "<0.001", formatC(p, format = "f", digits = 3)))
}

#' @export
print.tte_analysis_set <- function(x, ...) {
  cat("<tte_analysis_set> ", x$key, ": n = ", x$n, ", clusters = ",
      x$n_clusters, ", cases = ", x$cases, ", competing = ", x$competing,
      "\n", sep = "")
  invisible(x)
}
#' @export
print.tte_comparison <- function(x, ...) {
  cat("<tte_comparison> tolerance ", format(x$tol), "\n", sep = "")
  for (l in x$layers) {
    cat("  [", if (isTRUE(l$pass)) "PASS" else "FAIL", "] ", l$name,
        " - ", l$detail, "\n", sep = "")
  }
  invisible(x)
}
#' @export
print.tte_fit <- function(x, ...) {
  if (isTRUE(x$ok)) {
    cat("<tte_fit> ", x$key, " [", x$model, "] HR ",
        format_ci(x$est, x$lower, x$upper), ", p = ",
        signif(x$p, 3), "\n", sep = "")
  } else {
    cat("<tte_fit> ", x$key, " NOT FITTED: ", x$reason, "\n", sep = "")
  }
  invisible(x)
}
#' @export
print.tte_matched_set <- function(x, ...) {
  d <- x$data
  cat("<tte_matched_set>\n")
  cat("  pairs        : ", length(unique(d[["subclass"]])), "\n", sep = "")
  cat("  rows         : ", nrow(d), "\n", sep = "")
  cat("  prefiltered  : ", x$n_prefiltered, " participants entered the PS model\n",
      sep = "")
  cat("  cache key    : ", if (is.na(x$cache_key)) "<uncached>" else x$cache_key,
      "\n", sep = "")
  invisible(x)
}
#' @export
print.tte_outcome_registry <- function(x, ...) {
  d <- as.data.frame(x)[, c("key", "label", "model", "in_plot"), drop = FALSE]
  cat("<tte_outcome_registry> ", nrow(d), " outcome(s), ",
      sum(d$in_plot), " in plot\n", sep = "")
  print(d, row.names = FALSE)
  invisible(x)
}
#' @export
print.tte_params <- function(x, ...) {
  cat("<tte_params>\n")
  cat("  id / exposure : ", x$id_col, " / ", x$exposure_col, "\n", sep = "")
  cat("  time zero     : ", x$time_zero_col, "\n", sep = "")
  cat("  ps formula    : ", paste(deparse(x$ps_formula), collapse = " "),
      "\n", sep = "")
  cat("  calipers      : ps ", x$ps_caliper,
      if (length(x$match_calipers)) {
        paste0("; ", paste(names(x$match_calipers), x$match_calipers,
                           sep = " ", collapse = "; "))
      } else "", "\n", sep = "")
  cat("  year_days     : ", x$year_days,
      if (x$year_days == 365) "" else "  (baseline used 365)", "\n", sep = "")
  cat("  min_events    : ", x$min_events,
      if (x$min_events == 0) "  (guard inert)" else "", "\n", sep = "")
  cat("  followup_end  : ",
      if (is.null(x$followup_end)) "NULL (encoded upstream)"
      else if (inherits(x$followup_end, "Date")) format(x$followup_end)
      else paste0("column \"", x$followup_end, "\""), "\n", sep = "")
  invisible(x)
}
#' @export
print.tte_results <- function(x, ...) {
  ## `[.data.frame` keeps the tte_results class, so this method also receives
  ## column subsets like res[, c("key", "est")]. Fall back to printing whatever
  ## columns are present rather than assuming the full results schema.
  full <- c("label", "n", "n_clusters", "cases", "competing", "ci", "p", "q")
  cols <- intersect(full, names(x))
  if (!length(cols) || is.null(x$ok)) cols <- names(x)
  print(as.data.frame(x)[, cols, drop = FALSE], row.names = FALSE)
  if (!is.null(x$ok)) {
    nf <- sum(!x$ok)
    if (nf) cat("\n", nf, " outcome(s) not fitted; see the `reason` column.\n",
                sep = "")
  }
  invisible(x)
}
#' @export
print.tte_scenarios <- function(x, ...) {
  cat("<tte_scenarios> ", length(x$label), " arm(s), reference \"",
      x$reference, "\"\n", sep = "")
  for (i in seq_along(x$label)) {
    cat("  ", x$label[i], " [", x$strategy[[i]]$label, "] <- ",
        basename(x$source[i]), "\n", sep = "")
  }
  invisible(x)
}
#' @export
print.tte_strategy <- function(x, ...) {
  cat("<tte_strategy> ", x$label, " - ", x$mechanism, "\n", sep = "")
  invisible(x)
}
#' @export
summary.tte_matched_set <- function(object, ...) {
  cat("Available partners per exposed participant:\n")
  print(summary(object$match_num))
  invisible(object)
}
