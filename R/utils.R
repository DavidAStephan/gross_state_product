# Shared constants and small helpers ------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

STATES <- c("NSW", "VIC", "QLD", "SA", "WA", "TAS", "NT", "ACT")
STATE_NAMES <- c(
  NSW = "New South Wales", VIC = "Victoria", QLD = "Queensland",
  SA = "South Australia", WA = "Western Australia", TAS = "Tasmania",
  NT = "Northern Territory", ACT = "Australian Capital Territory"
)
BASES <- c("current", "cvm")
TOL <- 1e-6 # relative tolerance for the accounting constraints

DIR_RAW <- "data/raw"
DIR_PROC <- "data/processed"
DIR_OUT <- "output"

# ABS quarterly dates are stamped with the first day of the quarter's last
# month (2026-06-01 = June quarter 2026). Financial year = year ending June.
fy_of <- function(date) {
  as.integer(format(date, "%Y")) + as.integer(as.integer(format(date, "%m")) > 6)
}
# Quarter within the financial year: Sep = 1, Dec = 2, Mar = 3, Jun = 4
fyq_of <- function(date) {
  m <- as.integer(format(date, "%m"))
  c(3L, 3L, 3L, 4L, 4L, 4L, 1L, 1L, 1L, 2L, 2L, 2L)[m]
}
# Quarter-end-month date for a monthly date
qdate_of <- function(date) {
  m <- as.integer(format(date, "%m"))
  y <- as.integer(format(date, "%Y"))
  as.Date(sprintf("%d-%02d-01", y, ceiling(m / 3) * 3))
}

# Stop with an instruction the user can act on. Never fall back to fake data.
stop_data <- function(what, file_hint) {
  stop(
    "\nDATA PROBLEM: ", what,
    "\nPlease supply: ", file_hint,
    "\nNo placeholder data will be substituted.", call. = FALSE
  )
}

# Wide matrix (rows = dates, cols = STATES) from a long data frame
to_wide <- function(df, value, dates = NULL) {
  w <- df |>
    select(date, state, value = all_of(value)) |>
    pivot_wider(names_from = state, values_from = value) |>
    arrange(date)
  if (!is.null(dates)) w <- w[match(dates, w$date), ]
  m <- as.matrix(w[, STATES])
  rownames(m) <- as.character(w$date)
  m
}

rel_err <- function(x, target) abs(x - target) / pmax(abs(target), 1)

pct <- function(x) 100 * (x / dplyr::lag(x) - 1)

msg <- function(...) cat(format(Sys.time(), "%H:%M:%S"), "|", ..., "\n")

save_proc <- function(x, name) {
  dir.create(DIR_PROC, showWarnings = FALSE, recursive = TRUE)
  saveRDS(x, file.path(DIR_PROC, paste0(name, ".rds")))
}
read_proc <- function(name) {
  f <- file.path(DIR_PROC, paste0(name, ".rds"))
  if (!file.exists(f)) stop("Missing ", f, " - run the earlier numbered scripts first.", call. = FALSE)
  readRDS(f)
}
