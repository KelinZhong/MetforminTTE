# ---------------------------------------------------------------------------
# Configuration: the four objects you build before running anything.
#
# Nothing about the study is hardcoded. Column names, the propensity formula,
# calipers and the follow-up rule live in tte_params(); the outcome list lives
# in tte_outcomes(); the estimand lives in tte_strategy(); the arms live in
# tte_scenarios().
# ---------------------------------------------------------------------------

`%||%` <- function(a, b) if (is.null(a)) b else a

#' Analysis parameters
#'
#' Every column name and rule the pipeline needs. `id_col`, `exposure_col`,
#' `ps_formula` and `time_zero_col` are required with no defaults, so a study
#' specific column name can never be inherited by accident.
#'
#' @param id_col Participant identifier column.
#' @param exposure_col Binary exposure column (1 = new user, 0 = eligible
#'   control).
#' @param ps_formula Propensity-score formula; its left-hand side must be
#'   `exposure_col`.
#' @param time_zero_col Date column supplying time zero. Both members of a
#'   matched pair take the **exposed** member's value.
#' @param death_status_col,death_date_col Mortality indicator and date.
#'   Required for `model = "crr"` outcomes, where death is the competing
#'   event.
#' @param match_calipers Named numeric vector of exact-scale calipers applied
#'   jointly with the propensity caliper, e.g. `c(hba1c = 0.25)`. May be empty.
#' @param ps_caliper Propensity-score caliper. Default 0.05.
#' @param year_days Days per year. Default 365 — **not** 365.25, matching the
#'   original analysis.
#' @param min_events Minimum events before an outcome is fitted. Default 0,
#'   which disables the guard. Raise it only after equivalence with a baseline
#'   is established: a non-zero value suppresses fits the original produced.
#' @param followup_end End of administrative follow-up: a `Date`, a column
#'   name, or `NULL` meaning each outcome's date column already carries the
#'   event-or-censoring date for every row.
#' @param prev_affected_level Level of the `*_prev` columns denoting presence
#'   of the condition, e.g. `"1"` for 0/1 coding or `"2"` for 1/2 coding.
#'   `NULL` (default) leaves the coding untouched. Read this off
#'   [tte_validate()]'s tables — guessing it inverts every exclusion in the
#'   study without producing a single error.
#' @param recode Function applied to the cohort after reading, or the string
#'   `"ukb"` for this study's UK Biobank category collapsing. Default
#'   `identity`.
#' @param reader Optional function overriding file-format dispatch, e.g.
#'   `arrow::read_parquet`.
#'
#' @return An object of class `tte_params`.
#' @examples
#' params <- tte_params(
#'   id_col = "eid", exposure_col = "met_user",
#'   ps_formula = met_user ~ age + sex + bmi,
#'   time_zero_col = "assessment_date",
#'   death_status_col = "death_status", death_date_col = "date_of_death",
#'   match_calipers = c(hba1c = 0.25),
#'   ps_caliper = 0.05)
#' params
#' @export
tte_params <- function(id_col,
                       exposure_col,
                       ps_formula,
                       time_zero_col,
                       death_status_col = NULL,
                       death_date_col = NULL,
                       match_calipers = numeric(0),
                       ps_caliper = 0.05,
                       year_days = 365,
                       min_events = 0,
                       followup_end = NULL,
                       prev_affected_level = NULL,
                       recode = identity,
                       reader = NULL) {

  if (missing(id_col) || !is.character(id_col) || length(id_col) != 1L) {
    stop("`id_col` is required and must be a single column name.", call. = FALSE)
  }
  if (missing(exposure_col) || !is.character(exposure_col) ||
      length(exposure_col) != 1L) {
    stop("`exposure_col` is required and must be a single column name.",
         call. = FALSE)
  }
  if (missing(ps_formula) || !inherits(ps_formula, "formula")) {
    stop("`ps_formula` is required and must be a formula.", call. = FALSE)
  }
  if (missing(time_zero_col) || !is.character(time_zero_col) ||
      length(time_zero_col) != 1L) {
    stop("`time_zero_col` is required and must be a single column name.",
         call. = FALSE)
  }

  lhs <- all.vars(ps_formula[[2L]])
  if (length(lhs) != 1L || !identical(lhs, exposure_col)) {
    stop("`ps_formula` left-hand side must be `exposure_col` (\"",
         exposure_col, "\"), not \"", paste(lhs, collapse = " + "), "\".",
         call. = FALSE)
  }

  if (length(match_calipers) &&
      (is.null(names(match_calipers)) || any(!nzchar(names(match_calipers))))) {
    stop("`match_calipers` must be a *named* numeric vector, e.g. ",
         "c(hba1c = 0.25).", call. = FALSE)
  }
  if (!is.numeric(ps_caliper) || length(ps_caliper) != 1L || ps_caliper <= 0) {
    stop("`ps_caliper` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(year_days) || length(year_days) != 1L || year_days <= 0) {
    stop("`year_days` must be a single positive number.", call. = FALSE)
  }
  if (!is.numeric(min_events) || length(min_events) != 1L || min_events < 0) {
    stop("`min_events` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.null(followup_end) &&
      !inherits(followup_end, "Date") &&
      !(is.character(followup_end) && length(followup_end) == 1L)) {
    stop("`followup_end` must be a Date, a single column name, or NULL.",
         call. = FALSE)
  }
  if (is.character(recode)) {
    recode <- switch(recode,
      ukb = tte_recode_ukb,
      identity = identity,
      stop("`recode` must be a function, \"ukb\", or \"identity\".",
           call. = FALSE))
  }
  if (!is.function(recode)) {
    stop("`recode` must be a function, \"ukb\", or \"identity\".", call. = FALSE)
  }

  structure(
    list(id_col = id_col,
         exposure_col = exposure_col,
         ps_formula = ps_formula,
         time_zero_col = time_zero_col,
         death_status_col = death_status_col,
         death_date_col = death_date_col,
         match_calipers = match_calipers,
         ps_caliper = ps_caliper,
         year_days = year_days,
         min_events = min_events,
         followup_end = followup_end,
         prev_affected_level = prev_affected_level,
         recode = recode,
         reader = reader),
    class = "tte_params"
  )
}

#' Outcome registry
#'
#' One row per outcome. Column names are supplied per outcome, so irregular
#' naming is configuration rather than a code special case — an outcome whose
#' date column is `any_can_date` instead of `any_can_diag_date` needs no code
#' change.
#'
#' Elementary outcomes are vectorised: pass vectors and scalars recycle.
#' Supplying `components` instead builds a single **derived** (composite)
#' outcome from elementary keys.
#'
#' `from` records, per component, whether the composite consumes that
#' component's competing-risk-recoded indicator (`"cr"`) or its raw incidence
#' indicator (`"raw"`). This exists so a known inconsistency in an original
#' analysis can be reproduced as recorded configuration and corrected later
#' with a one-argument edit rather than a code change.
#'
#' @param key Short identifier; names the derived columns `time_to_<key>` and
#'   `<key>_incid_cr`, and merges arms in combined forest plots. Must be
#'   unique.
#' @param label Display label.
#' @param prev Baseline-prevalence column used to exclude pairs. For a derived
#'   outcome, `NA` means "union of the components' prevalence flags".
#' @param status Baseline-status column. Used only by censoring; outcome
#'   derivation never sees it.
#' @param incid First-occurrence indicator column.
#' @param date Diagnosis-date column.
#' @param model `"crr"` (Fine-Gray, death competing) or `"cox"` (mortality
#'   itself, which has no competing event).
#' @param in_plot Whether the outcome appears in forest plots. `FALSE` keeps it
#'   in the results table but out of the plot, replacing positional row drops.
#' @param components Character vector of elementary keys. Supplying this makes
#'   a derived outcome; `key` must then be length 1.
#' @param from `"cr"` or `"raw"`, length 1 (recycled) or `length(components)`.
#' @param dates Date columns for time-to-first-event in a derived outcome.
#'   Defaults to the components' own date columns.
#'
#' @return A `tte_outcome_registry` (a `data.table`). Combine pieces with
#'   `rbind()`.
#' @examples
#' keys <- c("dementia", "stroke", "diabetes")
#' elementary <- tte_outcomes(
#'   key    = keys,
#'   label  = c("Dementia", "Stroke", "Type 2 diabetes"),
#'   prev   = paste0(keys, "_prev"),
#'   status = paste0(keys, "_status"),
#'   incid  = paste0(keys, "_incid"),
#'   date   = paste0(keys, "_diag_date"),
#'   in_plot = keys != "diabetes")
#'
#' ## irregular column names are absorbed here, not special-cased in code
#' cancer <- tte_outcomes("any_can", "Any cancer", "any_can_prev",
#'                        "any_can_status", "any_can_incid",
#'                        date = "any_can_date")
#'
#' ## mortality has no competing event, so it is a Cox model
#' mortality <- tte_outcomes("mortality", "All-cause mortality",
#'                           status = "death_status", incid = "death_status",
#'                           date = "date_of_death", model = "cox")
#'
#' ## a composite; prev defaults to the union of its components' flags
#' composite <- tte_outcomes("composite", "Composite age-related disease",
#'                           components = c("dementia", "stroke"), from = "cr")
#'
#' outcomes <- rbind(elementary, cancer, mortality, composite)
#' outcomes
#' @export
tte_outcomes <- function(key, label, prev = NA_character_,
                         status = NA_character_, incid = NA_character_,
                         date = NA_character_, model = "crr", in_plot = TRUE,
                         components = NULL, from = "cr", dates = NULL) {

  if (!is.null(components)) {
    return(.derived_outcome(key, label, components, from, dates, prev,
                            in_plot))
  }

  n <- length(key)
  if (!n) stop("`key` must have at least one element.", call. = FALSE)
  rec <- function(x, nm) {
    if (length(x) == 1L) x <- rep(x, n)
    if (length(x) != n) {
      stop("`", nm, "` must be length 1 or length(key) = ", n, ".",
           call. = FALSE)
    }
    x
  }
  label  <- rec(label,  "label")
  prev   <- rec(prev,   "prev")
  status <- rec(status, "status")
  incid  <- rec(incid,  "incid")
  date   <- rec(date,   "date")
  model  <- rec(model,  "model")
  in_plot <- rec(in_plot, "in_plot")

  bad <- setdiff(unique(model), c("crr", "cox"))
  if (length(bad)) {
    stop("`model` must be \"crr\" or \"cox\"; got: ",
         paste(bad, collapse = ", "), call. = FALSE)
  }
  if (anyDuplicated(key)) {
    stop("`key` values must be unique; duplicated: ",
         paste(unique(key[duplicated(key)]), collapse = ", "), call. = FALSE)
  }

  ## Built via as.data.table(list(...)): calling data.table() directly would
  ## swallow our `key` column as its own reserved `key` argument (which sets
  ## index columns) and fail. The three list columns are created here too:
  ## assigning a list of NULLs to a data.table with `$<-` is interpreted as
  ## column REMOVAL, so they must exist from construction.
  out <- data.table::as.data.table(list(
    key = as.character(key),
    label = as.character(label),
    prev = as.character(prev),
    status = as.character(status),
    incid = as.character(incid),
    date = as.character(date),
    model = as.character(model),
    in_plot = as.logical(in_plot),
    components = vector("list", n),
    from = vector("list", n),
    dates = vector("list", n)
  ))

  class(out) <- unique(c("tte_outcome_registry", class(out)))
  out
}
#' @noRd
.derived_outcome <- function(key, label, components, from, dates, prev,
                             in_plot) {

  if (length(key) != 1L) {
    stop("a derived outcome must be built one at a time (length(key) == 1).", call. = FALSE)
  }
  if (!length(components)) {
    stop("`components` must name at least one elementary outcome.",
         call. = FALSE)
  }
  if (length(from) == 1L) from <- rep(from, length(components))
  if (length(from) != length(components)) {
    stop("`from` must be length 1 or length(components).", call. = FALSE)
  }
  bad <- setdiff(unique(from), c("cr", "raw"))
  if (length(bad)) {
    stop("`from` must be \"cr\" or \"raw\"; got: ", paste(bad, collapse = ", "),
         call. = FALSE)
  }

  ## as.data.table(list(...)) for the same reason as in tte_outcomes(): the
  ## `key` column name collides with data.table()'s reserved `key` argument.
  out <- data.table::as.data.table(list(
    key = as.character(key),
    label = as.character(label),
    prev = if (is.null(prev) || all(is.na(prev))) NA_character_ else as.character(prev)[1],
    status = NA_character_,
    incid = NA_character_,
    date = NA_character_,
    model = "derived",
    in_plot = as.logical(in_plot),
    components = list(as.character(components)),
    from = list(as.character(from)),
    dates = list(if (is.null(dates)) NULL else as.character(dates))
  ))

  class(out) <- unique(c("tte_outcome_registry", class(out)))
  out
}

#' Analysis strategy
#'
#' A strategy is the estimand plus the mechanism that produces it. Three
#' concepts are kept separate on purpose: the *label* is what you claim
#' (intention-to-treat, per-protocol), the *mechanism* is whatever rewrites
#' event times, and the *data* driving it is an argument. Collapsing these is
#' how "per-protocol" comes to mean whichever censoring rule someone wrote
#' last.
#'
#' @section Processor contract:
#' ```
#' in:  surv - time zero set, outcomes NOT yet derived
#' out: surv - status / incidence / date columns rewritten as needed
#'      may add a `censor_date` column; follow-up is capped at it
#'      output MUST NOT contain any time_to_* or *_incid_cr column
#' ```
#' The last line is enforced; derivation runs afterwards. A processor that
#' removes events must also write `censor_date`, or follow-up will not shorten
#' to match and the participant gains person-time they never contributed.
#'
#' @section Per-protocol assumptions:
#' Controls inherit their matched exposed participant's censoring date, so
#' pair follow-up stays symmetric — but a control's censoring is driven by
#' someone else's prescribing history. Administrative censoring at
#' discontinuation without weighting identifies the per-protocol effect only
#' under non-informative discontinuation; if discontinuation is informative,
#' take [tte_analysis_data()] and fit weighted models yourself.
#'
#' @param type `"itt"`, `"pp"`, or `"custom"`. A function may be passed
#'   directly as shorthand for a custom strategy.
#' @param rx Dispensing records (path or data.frame) for `"pp"`. Read eagerly,
#'   so malformed dates fail here rather than mid-run.
#' @param grace_days Days added to the last dispensing before censoring.
#' @param end_date Optional cap (administrative end of dispensing data).
#' @param label Strategy label; used in filenames and plot legends. Defaults
#'   to `"ITT"`, `"PP"` or `"custom"`.
#' @param process For `"custom"`, a function of `(surv, matched, params)`.
#' @param missing_rx Behaviour for an exposed participant with no dispensing
#'   records: `"uncensored"` (default, believed to reproduce the original),
#'   `"censor_at_zero"`, or `"error"`. Confirm against your baseline.
#' @param eid_col,date_col Column names in the dispensing data.
#'
#' @return An object of class `tte_strategy`.
#' @examples
#' itt <- tte_strategy("itt")
#' itt
#'
#' rx <- data.frame(eid = c("P001", "P001", "P002"),
#'                  issue_date = as.Date(c("2011-01-01", "2012-06-01",
#'                                         "2011-03-01")))
#' pp <- tte_strategy("pp", rx = rx, grace_days = 180)
#' pp
#'
#' ## a marker strategy proves the hook runs before derivation, not merely
#' ## that it exists
#' marker <- tte_strategy(function(surv, matched, params) {
#'   surv$oa_date <- surv$time_zero + 365
#'   surv
#' })
#' marker$label
#' @export
tte_strategy <- function(type = c("itt", "pp", "custom"), rx = NULL,
                         grace_days = 180, end_date = NULL, label = NULL,
                         process = NULL,
                         missing_rx = c("uncensored", "censor_at_zero", "error"),
                         eid_col = "eid", date_col = "issue_date") {

  if (is.function(type)) { process <- type; type <- "custom" }
  type <- match.arg(type)
  missing_rx <- match.arg(missing_rx)

  if (identical(type, "custom")) {
    if (!is.function(process)) {
      stop("A custom strategy needs `process`, a function of ",
           "(surv, matched, params).", call. = FALSE)
    }
    if (length(formals(process)) < 3L) {
      stop("`process` must accept three arguments: (surv, matched, params).",
           call. = FALSE)
    }
    return(structure(list(label = label %||% "custom", process = process,
                          mechanism = "custom"), class = "tte_strategy"))
  }

  if (identical(type, "itt")) {
    return(structure(list(label = label %||% "ITT",
                          process = function(surv, matched, params) surv,
                          mechanism = "no censoring"),
                     class = "tte_strategy"))
  }

  ## type == "pp"
  if (is.null(rx)) {
    stop("A per-protocol strategy needs `rx`, the dispensing records.",
         call. = FALSE)
  }
  if (!is.null(end_date) && !inherits(end_date, "Date")) {
    stop("`end_date` must be a Date or NULL.", call. = FALSE)
  }

  d <- tte_read(rx)
  miss <- setdiff(c(eid_col, date_col), names(d))
  if (length(miss)) {
    stop("Dispensing file missing column(s): ", paste(miss, collapse = ", "),
         call. = FALSE)
  }

  parsed <- .as_date(d[[date_col]])
  bad <- is.na(parsed) & !is.na(d[[date_col]])
  if (any(bad)) {
    stop(sum(bad), " unparseable dispensing date(s), e.g. ",
         paste(utils::head(unique(d[[date_col]][bad]), 5), collapse = ", "),
         ". Fix the file or pass a pre-parsed data.frame.", call. = FALSE)
  }

  rx_dt <- data.table::data.table(eid = d[[eid_col]], issue_date = parsed)
  n_before <- nrow(rx_dt)
  rx_dt <- unique(rx_dt)
  if (nrow(rx_dt) < n_before) {
    message("tte_strategy(): removed ", n_before - nrow(rx_dt),
            " duplicate dispensing rows.")
  }
  rx_dt <- rx_dt[order(rx_dt$eid, rx_dt$issue_date), ]

  structure(list(
    label = label %||% "PP",
    process = function(surv, matched, params) {
      .pp_process(surv, matched, params, rx_dt = rx_dt,
                  grace_days = grace_days, end_date = end_date,
                  missing_rx = missing_rx)
    },
    mechanism = paste0("censors at last dispensing + ", grace_days, "d",
                       if (!is.null(end_date))
                         paste0(", capped at ", format(end_date)) else ""),
    rx = rx_dt, grace_days = grace_days, end_date = end_date,
    missing_rx = missing_rx), class = "tte_strategy")
}

#' Scenario (arm) table
#'
#' @param sources **Named** character vector of cohort file paths; names are
#'   arm labels used for output filenames and plot legends. Any number of arms.
#' @param strategy List of [tte_strategy()] objects, one per source.
#' @param reference Label of the arm whose ordering governs combined forest
#'   plots. Defaults to the first arm — never a hardcoded literal.
#'
#' @details Arms sharing a source path (and matching-relevant parameters)
#'   share one cached matched set, so an intention-to-treat arm and a
#'   per-protocol arm built from the same file have *identical* pairs by
#'   construction rather than by coincidence.
#'
#' @return An object of class `tte_scenarios`.
#' @examples
#' rx <- data.frame(eid = "P001", issue_date = as.Date("2011-01-01"))
#' scen <- tte_scenarios(
#'   sources  = c(`1yr ITT` = "cohort_1yr.csv", `1yr PP` = "cohort_1yr.csv"),
#'   strategy = list(tte_strategy("itt"), tte_strategy("pp", rx = rx)),
#'   reference = "1yr ITT")
#' scen
#' @export
tte_scenarios <- function(sources, strategy, reference = names(sources)[1L]) {

  if (is.null(names(sources)) || any(!nzchar(names(sources)))) {
    stop("`sources` must be a *named* character vector; names are arm labels.",
         call. = FALSE)
  }
  if (anyDuplicated(names(sources))) {
    stop("Arm labels must be unique.", call. = FALSE)
  }
  if (!is.list(strategy)) strategy <- list(strategy)
  if (length(strategy) != length(sources)) {
    stop("`strategy` must have one entry per source (", length(sources), ").",
         call. = FALSE)
  }
  ok <- vapply(strategy, inherits, logical(1), what = "tte_strategy")
  if (!all(ok)) {
    stop("Every `strategy` entry must be a tte_strategy object; bad entries: ",
         paste(which(!ok), collapse = ", "), call. = FALSE)
  }
  if (!reference %in% names(sources)) {
    stop("`reference` must be one of the arm labels: ",
         paste(names(sources), collapse = ", "), call. = FALSE)
  }

  structure(
    list(label = names(sources),
         source = unname(as.character(sources)),
         strategy = strategy,
         reference = reference),
    class = "tte_scenarios"
  )
}

#' @noRd
.pp_process <- function(surv, matched, params, rx_dt, grace_days, end_date,
                        missing_rx) {

  id_col <- params$id_col
  exp_col <- params$exposure_col
  outcomes <- attr(surv, "outcomes")

  ids <- surv[[id_col]]
  exposed <- !is.na(surv[[exp_col]]) & surv[[exp_col]] == 1

  ## per-participant dispensing summary (descriptives, and the natural
  ## covariate set for anyone later modelling discontinuation)
  agg <- rx_dt[rx_dt$eid %in% ids[exposed], ]
  summ <- NULL
  if (nrow(agg)) {
    sp <- split(agg$issue_date, agg$eid)
    summ <- data.frame(
      id = names(sp),
      n_rx = vapply(sp, length, integer(1)),
      first_rx = as.Date(vapply(sp, function(x) as.numeric(min(x)), numeric(1)),
                         origin = "1970-01-01"),
      last_rx = as.Date(vapply(sp, function(x) as.numeric(max(x)), numeric(1)),
                        origin = "1970-01-01"),
      max_gap_days = vapply(sp, function(x) {
        if (length(x) < 2L) 0 else max(as.numeric(diff(sort(x))))
      }, numeric(1)),
      stringsAsFactors = FALSE)
    summ$duration_days <- as.numeric(summ$last_rx - summ$first_rx)
    rownames(summ) <- NULL
  }

  last_rx <- rep(as.Date(NA), nrow(surv))
  if (!is.null(summ)) {
    hit <- match(as.character(ids), summ$id)
    last_rx <- summ$last_rx[hit]
  }

  ## exposed participants absent from the dispensing file
  zero <- exposed & is.na(last_rx)
  if (any(zero)) {
    if (identical(missing_rx, "error")) {
      stop(sum(zero), " exposed participant(s) have no dispensing records. ",
           "Set `missing_rx` to \"uncensored\" or \"censor_at_zero\".",
           call. = FALSE)
    }
    message("tte_strategy(): ", sum(zero),
            " exposed participant(s) with no dispensing records -> ",
            missing_rx, ".")
  }

  censor <- last_rx + grace_days
  if (!is.null(end_date)) censor[!is.na(censor)] <-
    pmin(censor[!is.na(censor)], end_date)
  if (identical(missing_rx, "censor_at_zero")) {
    censor[zero] <- .as_date(surv[["time_zero"]])[zero]
  }

  ## propagate the exposed member's date to their matched control
  sc <- surv[["subclass"]]
  user_date <- stats::setNames(censor[exposed], as.character(sc[exposed]))
  censor <- unname(user_date[as.character(sc)])
  surv[["censor_date"]] <- censor

  ## rewrite status / incidence / date for every elementary outcome
  elem <- outcomes[outcomes$model != "derived", ]
  for (i in seq_len(nrow(elem))) {
    surv <- apply_censoring(surv,
                            key = elem$key[i],
                            status = elem$status[i],
                            incid = elem$incid[i],
                            date = elem$date[i],
                            censor_date = censor)
  }
  if (!is.null(params$death_status_col) && !is.null(params$death_date_col)) {
    surv <- apply_censoring(surv, key = "death",
                            status = params$death_status_col,
                            incid = params$death_status_col,
                            date = params$death_date_col,
                            censor_date = censor)
  }

  attr(surv, "rx_summary") <- summ
  surv
}
