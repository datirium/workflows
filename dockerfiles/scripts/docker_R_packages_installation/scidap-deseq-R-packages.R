# install_packages.R

# List of CRAN packages to install, including dependencies
cran_packages <- c(
  "googledrive", "googlesheets4", "httr", "ragg", "rvest", "xml2",
  "locfit", "tidyverse", "patchwork", "argparse", "data.table",
  "pheatmap", "ggrepel", "htmlwidgets", "devtools", "kableExtra",
  "systemfonts", "svglite"
)

# List of Bioconductor packages to install
bioc_packages <- c(
  "BiocParallel", "limma", "EnhancedVolcano", "hopach", "cmapR", "XVector",
  "S4Arrays", "SparseArray", "GenomicRanges", "DelayedArray", "DESeq2", "sva"
)

# Function to install and check CRAN packages
install_and_check_cran <- function(packages) {
  install.packages(packages, repos = "https://cloud.r-project.org")
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("Package", pkg, "failed to install."))
    }
  }
}

# Function to install and check Bioconductor packages
install_and_check_bioc <- function(packages) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = "https://cloud.r-project.org")
  }
  BiocManager::install(version = "3.19", update = FALSE, ask = FALSE)
  BiocManager::install(packages, update = TRUE, ask = FALSE)
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE)) {
      stop(paste("Package", pkg, "failed to install."))
    }
  }
}

# Install CRAN packages
install_and_check_cran(cran_packages)

# Install Bioconductor packages
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

cat("All packages installed successfully.\n")