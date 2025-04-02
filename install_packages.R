# Install pak if not already available
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

# Define packages with versions
pkgs <- c(
  "GEOmap@2.5.11",
  "geomapdata@2.0.2", 
  "ReacTran@1.4.3.1",
  "geosphere@1.5.20",
  "deSolve@1.40",
  "magrittr@2.0.3",
  "maptools@1.1.8",
  "ggplot2@3.5.1",
  "cowplot@1.1.3",
  "animation@2.7",
  "RColorBrewer@1.1.3",
  "fields@16.3",
  "maps@3.4.2.1",
  "scales@1.3.0"
)

# Install with version control
pak::pkg_install(pkgs)
