# run_all.R ---------------------------------------------------------------------
# Run the whole pipeline in order and render the report.
#   Rscript run_all.R                      # fresh ABS download
#   GSP_REFRESH=false Rscript run_all.R    # reuse spreadsheets already in data/raw
#   GSP_TARGET=prorata Rscript run_all.R   # literal per-year national adjustment
scripts <- c("01_data.R", "02_indicators.R", "03_baselines.R", "04_method_a.R", "04b_bottom_up_test.R", "05_method_b.R", "05b_method_c.R", "06_validate.R", "07_dashboard_export.R")
for (s in scripts) {
  cat("\n==================== ", s, " ====================\n", sep = "")
  status <- system2(file.path(R.home("bin"), "Rscript"), s)
  if (status != 0) stop("Stopped: ", s, " failed. Fix the problem above and rerun.", call. = FALSE)
}
cat("\n==================== report.qmd ====================\n")
# Render next to the source and then move: pointing Quarto's --output-dir at
# output/ would delete the CSVs already there.
quarto::quarto_render("report.qmd", quiet = TRUE)
dir.create("output", showWarnings = FALSE)
file.rename("report.html", file.path("output", "report.html"))
cat("Done. See output/quarterly_gsp.csv and output/report.html\n")
