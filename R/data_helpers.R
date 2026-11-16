# Data access helpers ----------------------------------------------------------
# All ABS data come through readabs. Every series used is selected by an
# explicit table + series-name pattern that must match exactly one series ID,
# and is logged to the data dictionary.

suppressPackageStartupMessages(library(readabs))

.dict <- new.env()
.dict$rows <- list()

dict_add <- function(concept, state, basis, freq, d) {
  .dict$rows[[length(.dict$rows) + 1]] <- tibble::tibble(
    concept = concept, state = state, price_basis = basis, frequency = freq,
    catalogue = d$cat[1], table = d$table_title[1], series = gsub(" +", " ", d$series[1]),
    series_type = d$series_type[1], series_id = d$series_id[1], unit = d$unit[1],
    first_obs = min(d$date[!is.na(d$value)]), last_obs = max(d$date[!is.na(d$value)])
  )
}
dict_get <- function() dplyr::bind_rows(.dict$rows)

# Download (or reuse local copies of) ABS time series tables
get_abs <- function(cat_no, tables, refresh = TRUE) {
  dir.create(DIR_RAW, showWarnings = FALSE, recursive = TRUE)
  out <- tryCatch(
    suppressMessages(read_abs(cat_no,
      tables = tables, path = DIR_RAW,
      check_local = !refresh, show_progress_bars = FALSE
    )),
    error = function(e) e
  )
  if (inherits(out, "error") || nrow(out) == 0) {
    stop_data(
      paste0("could not download ABS catalogue ", cat_no, " tables ",
             paste(tables, collapse = ", "), " (", conditionMessage(out), ")"),
      paste0("the ABS time series spreadsheets for catalogue ", cat_no, ", tables ",
             paste(tables, collapse = ", "), ", saved into ", DIR_RAW,
             "/ (then rerun with GSP_REFRESH=false)")
    )
  }
  out$cat <- cat_no
  out
}

# Pick exactly one series; fail loudly otherwise
pick <- function(df, table_rx, series_rx, type = "Seasonally Adjusted",
                 concept, state = NA_character_, basis = NA_character_, freq = "Q") {
  d <- df[grepl(table_rx, df$table_no) & grepl(series_rx, df$series) & df$series_type == type, ]
  ids <- unique(d$series_id)
  if (length(ids) != 1) {
    stop_data(
      sprintf("expected exactly one series for '%s' (%s, %s) but found %d [table /%s/, series /%s/, type %s]",
              concept, state, basis, length(ids), table_rx, series_rx, type),
      sprintf("the ABS table matching '%s' containing that series, or update the pattern in 01_data.R", table_rx)
    )
  }
  d <- d[!is.na(d$value), ]
  dict_add(concept, state, basis, freq, d)
  d[order(d$date), c("date", "value")]
}

# Seasonally adjust a monthly series with a robust STL (log scale when the
# series is strictly positive). Used only where the ABS publishes no SA series.
stl_sa <- function(date, value) {
  o <- order(date)
  date <- date[o]; value <- value[o]
  st <- c(as.integer(format(date[1], "%Y")), as.integer(format(date[1], "%m")))
  use_log <- all(value > 0)
  x <- stats::ts(if (use_log) log(value) else value, start = st, frequency = 12)
  fit <- stats::stl(x, s.window = 13, robust = TRUE)
  adj <- as.numeric(x - fit$time.series[, "seasonal"])
  tibble::tibble(date = date, value = if (use_log) exp(adj) else adj)
}

# Sum complete months to quarters (drops incomplete quarters)
m_to_q <- function(d) {
  d |>
    dplyr::mutate(qdate = qdate_of(date)) |>
    dplyr::group_by(date = qdate) |>
    dplyr::summarise(n = dplyr::n(), value = sum(value), .groups = "drop") |>
    dplyr::filter(n == 3) |>
    dplyr::select(date, value)
}

# Queensland Treasury quarterly state accounts (cross-check only) ----------------
# Publication is suspended and the tables have been removed from the live QGSO
# site. Sources tried, newest issue first:
#   1. September quarter 2020 issue (the last one) at its original QGSO address;
#   2. the Internet Archive copy of that same file;
#   3. March quarter 2019 issue on QGSO's own archive site (archive.qgso.qld.gov.au).
# A workbook already in data/raw/qld is used as-is (newest issue first).
QLD_SOURCES <- c(
  "qld-state-accounts-tables-202009.xlsx" = "https://www.qgso.qld.gov.au/issues/3371/qld-state-accounts-tables-202009.xlsx",
  "qld-state-accounts-tables-202009.xlsx" = "https://web.archive.org/web/2022id_/https://www.qgso.qld.gov.au/issues/3371/qld-state-accounts-tables-202009.xlsx",
  "qld-state-accounts-tables-201903.xlsx" = "https://archive.qgso.qld.gov.au/products/tables/qld-state-accounts-tables/qld-state-accounts-tables-201903.xlsx"
)
get_qld_workbook <- function() {
  dir <- file.path(DIR_RAW, "qld")
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ok <- function(f) file.exists(f) && file.size(f) > 20000 &&
    !inherits(try(readxl::excel_sheets(f), silent = TRUE), "try-error")
  for (f in unique(file.path(dir, names(QLD_SOURCES)))) if (ok(f)) return(f)
  for (i in seq_along(QLD_SOURCES)) {
    dest <- file.path(dir, names(QLD_SOURCES)[i])
    suppressWarnings(try(utils::download.file(QLD_SOURCES[[i]], dest, mode = "wb", quiet = TRUE,
                                              headers = c("User-Agent" = "Mozilla/5.0")), silent = TRUE))
    if (ok(dest)) return(dest)
    if (file.exists(dest)) unlink(dest)
  }
  stop_data(
    "could not download the Queensland Treasury quarterly state accounts workbook",
    paste0("'qld-state-accounts-tables-202009.xlsx' (Queensland State Accounts tables, September quarter 2020) or ",
           "'qld-state-accounts-tables-201903.xlsx' (March quarter 2019), saved in ", dir, "/")
  )
}

# Seasonally adjusted quarterly Queensland GSP ($m), chain volume and current prices
read_qld <- function(file) {
  sheets <- setdiff(readxl::excel_sheets(file), "Index")
  titles <- vapply(sheets, function(sh) {
    a1 <- suppressMessages(readxl::read_excel(file, sheet = sh, range = "A1", col_names = FALSE))
    if (ncol(a1) == 0 || nrow(a1) == 0) "" else as.character(a1[[1]])[1] # some issues contain an empty sheet
  }, character(1))
  titles[is.na(titles)] <- ""
  one <- function(rx, b) {
    sh <- sheets[grepl(rx, titles)]
    if (length(sh) != 1) stop_data(paste0("could not identify the sheet '", rx, "' in ", basename(file)),
                                   "the Queensland state accounts tables workbook in its published layout")
    x <- as.data.frame(suppressMessages(suppressWarnings(readxl::read_excel(file, sheet = sh, col_names = FALSE))))
    drow <- which(apply(x, 1, function(r) sum(suppressWarnings(as.numeric(r)) > 30000, na.rm = TRUE)) > 50)[1]
    grow <- which(trimws(x[[1]]) == "Gross state product")
    # the current-price table shows GSP twice (income and expenditure sides); they must agree
    if (length(grow) == 2) {
      g <- sapply(grow, function(i) suppressWarnings(as.numeric(unlist(x[i, ]))))
      if (max(abs(g[, 1] - g[, 2]), na.rm = TRUE) > 1) grow <- integer(0) else grow <- grow[1]
    }
    if (is.na(drow) || length(grow) != 1) stop_data(paste("unexpected layout in sheet", sh, "of", basename(file)),
                                                    "the Queensland state accounts tables workbook in its published layout")
    cols <- which(!is.na(suppressWarnings(as.numeric(unlist(x[drow, ])))))
    d <- as.Date(as.numeric(unlist(x[drow, cols])), origin = "1899-12-30")
    tibble::tibble(date = as.Date(format(d, "%Y-%m-01")), basis = b, qld_treasury = as.numeric(unlist(x[grow, cols])))
  }
  dplyr::bind_rows(
    one("Queensland domestic production account, Seasonally adjusted, chain volume measure \\(a\\), \\$m", "cvm"),
    one("Queensland domestic production account, Seasonally adjusted, current prices, \\$m", "current")
  )
}
