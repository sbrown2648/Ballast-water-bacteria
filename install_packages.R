# install_packages.R
# Run this script to install all required packages into the renv library.
# After running, execute renv::snapshot() to update renv.lock.

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv")
}

# Install BiocManager first (needed for Bioconductor packages)
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# CRAN packages
cran_packages <- c(
  "ggplot2", "vegan", "plyr", "dplyr", "scales", "grid", "reshape2",
  "picante", "tidyr", "viridis", "patchwork", "RColorBrewer", "plotly",
  "seecolor", "tibble", "tidyverse", "kableExtra", "knitr", "rmarkdown",
  "devtools", "readr", "magrittr", "cowplot", "ggforce", "ggplotify",
  "ggprism", "magick", "svglite", "yaml"
)

install.packages(cran_packages)

# Bioconductor packages
bioc_packages <- c(
  "phyloseq", "DESeq2", "ComplexHeatmap", "ALDEx2",
  "SummarizedExperiment", "Biobase", "BiocGenerics",
  "metagenomeSeq", "edgeR", "lefser", "limma", "KEGGREST",
  "MicrobiomeStat"
)

BiocManager::install(bioc_packages)

# GitHub packages
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

github_packages <- c(
  "jbisanz/qiime2R",
  "mikemc/speedyseq",
  "david-barnett/microViz",
  "gmteunisse/fantaxtic",
  "gmteunisse/ggnested",
  "mahendra-mariadassou/phyloseq-extended",
  "pmartinezarbizu/pairwiseAdonis/pairwiseAdonis",
  "cafferychen777/ggpicrust2",
  "coolbutuseless/ggpattern"
)

for (pkg in github_packages) {
  remotes::install_github(pkg, upgrade = "never")
}

# ggh4x from CRAN
install.packages("ggh4x")

# After installation, snapshot the library
renv::snapshot(prompt = FALSE)

message("All packages installed. renv.lock has been updated.")
