# install_packages.R

# List of CRAN packages to install, including dependencies
cran_packages <- c(
  "googledrive", "googlesheets4", "httr", "ragg", "rvest", "xml2",
  "locfit", "tidyverse", "patchwork", "argparse", "data.table",
  "pheatmap", "ggrepel", "htmlwidgets", "devtools", "kableExtra",
  "systemfonts", "svglite", "gridExtra", "rlang", "stringr", 
  "RColorBrewer", "glue", "pryr", "plotly", "ggplot2", "future",
  "jsonlite", "magrittr", "reshape2", "conflicted", "logger", "params"
)

# List of Bioconductor packages to install
bioc_packages <- c(
  "BiocParallel", "limma", "EnhancedVolcano", "hopach", "cmapR", "XVector",
  "S4Arrays", "SparseArray", "GenomicRanges", "DelayedArray", "DESeq2", "sva"
)

# Function to install and check CRAN packages
install_and_check_cran <- function(packages) {
  install.packages(packages, repos = "https://cloud.r-project.org", dependencies = TRUE)
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      cat(paste("Warning: Package", pkg, "failed to install or load. Trying alternative method...\n"))
      try({
        install.packages(pkg, repos = "https://cloud.r-project.org", dependencies = TRUE)
        if (!require(pkg, character.only = TRUE)) {
          stop(paste("Package", pkg, "failed to install after retry."))
        }
      })
    } else {
      cat(paste("Package", pkg, "successfully installed and loaded.\n"))
    }
  }
}

# Function to install and check Bioconductor packages
install_and_check_bioc <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(version = "3.19", update = FALSE, ask = FALSE)
  BiocManager::install(packages, update = FALSE, ask = FALSE)
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      cat(paste("Warning: Bioconductor package", pkg, "failed to install or load. Trying alternative method...\n"))
      try({
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
        if (!require(pkg, character.only = TRUE)) {
          stop(paste("Bioconductor package", pkg, "failed to install after retry."))
        }
      })
    } else {
      cat(paste("Bioconductor package", pkg, "successfully installed and loaded.\n"))
    }
  }
}

# Install CRAN packages
cat("Installing CRAN packages...\n")
install_and_check_cran(cran_packages)

# Install Bioconductor packages
cat("Installing Bioconductor packages...\n")
install_and_check_bioc(bioc_packages)

# Install additional packages from GitHub
install.packages("devtools", repos = "https://cloud.r-project.org")
if (!require("devtools", character.only = TRUE)) {
  stop("Package devtools failed to install.")
}
# This package is required for the GlimmaV2 package to horizontal placement of the drop-down menu
install.packages(
  "https://bioconductor.org/packages/3.17/bioc/src/contrib/Glimma_2.10.0.tar.gz",
  repos = NULL,
  type = "source"
)
if (!require("Glimma", character.only = TRUE)) {
  stop("Package GlimmaV2 failed to install from GitHub.")
}

# Set global R options for better memory handling in scripts
cat("Setting up global R options for memory management...\n")
options(future.globals.maxSize = 4000 * 1024^2)  # 4GB max for global data
options(expressions = 5000)  # Increase expression stack size
options(gc.aggressiveness = 0)  # Default GC behavior
options(stringsAsFactors = FALSE)  # Don't convert strings to factors by default

cat("All packages installed successfully.\n")
cat("Memory management options have been configured.\n")