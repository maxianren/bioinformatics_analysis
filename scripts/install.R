if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}
library(clusterProfiler)

if (!require("org.Mm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Mm.eg.db")
}
library(org.Mm.eg.db)

if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!require("ComplexHeatmap", quietly = TRUE)) {
  install.packages("ComplexHeatmap", dependencies = TRUE)
}