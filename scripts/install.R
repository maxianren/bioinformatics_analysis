if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}

if (!require("GO.db", quietly = TRUE)) {
  BiocManager::install("GO.db")
}
if (!require("clusterProfiler", quietly = TRUE)) {
  BiocManager::install("clusterProfiler")
}

if (!require("org.Mm.eg.db", quietly = TRUE)) {
  BiocManager::install("org.Mm.eg.db")
}

if (!require("tidyverse", quietly = TRUE)) {
  install.packages("tidyverse")
}

if (!require("ComplexHeatmap", quietly = TRUE)) {
  tryCatch({
    install.packages("ComplexHeatmap", dependencies = TRUE)
  }, error = function(e) {
    if (grepl("not available for this version of R", e$message)) {
      if (!require("devtools", quietly = TRUE)) {
        install.packages("devtools", dependencies = TRUE)
      }
      library(devtools)
      install_github("jockergoo/ComplexHeatmap", dependencies = TRUE)
    } else {
      stop(e) 
    }
  })
}

if (!require("ggplot", quietly = TRUE)) {
  install.packages("ggplot", dependencies = TRUE)
}