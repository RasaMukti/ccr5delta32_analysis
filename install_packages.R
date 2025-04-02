# Install pak if not already available
if (!requireNamespace("pak", quietly = TRUE)) {
  install.packages("pak")
}

# Define packages with versions
pkgs <- c(
  "geomapdata@2.0.2", 
  "ReacTran@1.4.3.1",
  "geosphere@1.5.20",
  "deSolve@1.40",
  "magrittr@2.0.3",
  "cowplot@1.1.3",
  "animation@2.7",
  "RColorBrewer@1.1.3",
  "fields@16.3"
)

# Install with version control
pak::pkg_install(pkgs)
