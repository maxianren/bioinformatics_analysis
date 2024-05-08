if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")
s
BiocManager::install("clusterProfiler")

BiocManager::install("org.Mm.eg.db")

install.packages("tidyverse")
