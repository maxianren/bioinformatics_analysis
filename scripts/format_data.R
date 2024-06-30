library(R6)
library(tidyr)
library(dplyr)
library(readr)
library(readxl)
library(tidyverse)
library(org.Mm.eg.db)
library(AnnotationDbi)

GeneDataProcessor <- R6Class("GeneDataProcessor",
                         public = list(
                           directory = NULL,
                           pattern = NULL,

                           initialize = function(directory, pattern) {
                             self$directory <- directory
                             self$pattern <- pattern
                           },
                           
                           get_file_paths = function() {
                             file_paths <- list.files(self$directory, self$pattern, full.names = TRUE)
                             return(file_paths)
                           },
                           parse_file_paths = function() {
                             file_name_full <- self$get_file_paths()
                             file_name <- basename(file_name_full)
                             file_name_clean <- gsub(" .*", "", file_name)
                             
                             file_name_clean_cy <- file_name_clean
                             
                             sample_id <- stringr::str_extract(file_name_clean_cy, "(?<=_)\\d+(?=_)")
                             file_name_clean_cy <- stringr::str_remove(file_name_clean_cy, "_\\d+")
                             genotype <- stringr::str_extract(file_name_clean_cy, "(?<=_)\\w+(?=-)")
                             file_name_clean_cy <- stringr::str_remove(file_name_clean_cy, "_\\w+")
                             gender <- stringr::str_extract(file_name_clean_cy, "(?<=-)\\w(?=-)")
                             file_name_clean_cy <- stringr::str_remove(file_name_clean_cy, "-\\w")
                             treatment <- stringr::str_extract(file_name_clean_cy, "(?<=-).+")
                             
                             df <- data.frame(
                               file_name_full = file_name_full,
                               file_name_clean = file_name_clean,
                               sample_id = sample_id,
                               genotype = genotype,
                               gender = gender,
                               treatment = treatment,
                               stringsAsFactors = FALSE
                             ) 
                             return(df)
                           }
                         )
)

GeneFeatureExtractor <- R6Class("GeneFeatureExtractor",
                                inherit = GeneDataProcessor,
                                
                                public = list(
                                  feature = NULL,
                                  
                                  initialize = function(directory, pattern, feature) {
                                    super$initialize(directory, pattern)
                                    self$feature <- feature
                                  },
                                  
                                  extract_feature_pivot = function(df) {
                                    file_paths <- self$get_file_paths()
                                    updated_df <- df
                                    
                                    for (file_path in file_paths) {
                                      data_df <- read_csv(file_path)
                                      file_name <- gsub(pattern = " .*", replacement = "", x = basename(file_path))
                                      file_name_column <- paste(file_name, self$feature, sep = "_")
                                      
                                      merged_df <- updated_df %>%
                                        left_join(data_df, by = c("Name" = "Name")) %>%
                                        dplyr::select(Name, self$feature) %>%
                                        dplyr::rename(!!file_name_column := self$feature)
                                      
                                      updated_df[[file_name_column]] <- merged_df[[file_name_column]]
                                      updated_df <- updated_df %>%
                                        filter(!if_any(ends_with(paste0("_", self$feature)), is.na))
                                    }
                                    return(updated_df)
                                  },
                                  
                                  extract_feature_unpivot = function(df_pivot) {
                                    df_unpivot <- df_pivot %>%
                                      pivot_longer(
                                        cols = ends_with(paste0("_", self$feature)),
                                        names_to = "Sample",
                                        values_to = self$feature
                                      ) %>%
                                      dplyr::mutate(
                                        Sample = str_remove(Sample, paste0("_", self$feature)),
                                        Number = str_extract(Sample, "(?<=_)\\d+(?=_?)"),
                                        Sample = str_remove(Sample, "_\\d+_"),
                                        Gene = str_extract(Sample, "^[^\\-]+"),
                                        Sample = str_remove(Sample, "^[^\\-]+\\-"),
                                        Gender = str_sub(Sample, 1, 1),
                                        Treatment = str_sub(Sample, 3)
                                      ) %>%
                                      dplyr::select(-Sample)
                                    return(df_unpivot)
                                  },
                                  
                                  get_feature_hierarchy = function(suffix = paste0("_", self$feature)) {

                                    df <- self$parse_file_paths()
                                    df$file_name_clean <- paste0(df$file_name_clean, suffix)
                                    df <- df %>%
                                      arrange(desc(treatment), desc(genotype), desc(gender))
                                    return(df)
                                  }
                                )
)

DifferentialExpressionDataProcessor <- R6Class("DifferentialExpressionDataProcessor",
                                  public = list(
                                    file_path = NULL,
                                    df_diff_expr = NULL,
                                    
                                    initialize = function(file_path) {
                                      self$file_path <- file_path
                                      self$df_diff_expr <- read.csv(file_path, check.names = FALSE)
                                    },
                
                                    #filter diff_expr dataframe, return a df
                                    filter_df = function(log_2_FC, p_value) {
                                      df_filtered <- self$df_diff_expr %>%
                                        filter(!is.na(`P-value`) & abs(`Log? fold change`) > log_2_FC & `P-value` < p_value) 
                                      return(df_filtered)
                                    },

                                    #filter diff_expr dataframe, return a gene list
                                    filter_gene = function(log_2_FC, p_value) {
                                      genes_filtered <- self$filter_df(log_2_FC, p_value) %>%
                                        dplyr::select('Name')%>%
                                        pull(Name)
                                      return(genes_filtered)
                                    },
                                    
                                    get_filtered_genes_entrez_id = function(log_2_FC, p_value, gene_db) {
                                      genes <- self$filter_gene(log_2_FC, p_value)
                                      
                                      hs <- gene_db # org.Mm.eg.db
                                      my.symbols <- genes
                                      df_with_enter_id <- AnnotationDbi::select(hs,
                                                                                keys = my.symbols,
                                                                                columns = c("ENTREZID", "SYMBOL"),
                                                                                keytype = "SYMBOL") %>%
                                        dplyr::filter(!is.na(ENTREZID))
                                      
                                      return(df_with_enter_id)
                                    },
                                    # get top genes by log2 fild change, return df
                                    get_top_n_genes = function(top_n, log_2_FC, p_value) {
                                      df_filtered <- self$filter_df(log_2_FC, p_value)%>%
                                        arrange(desc(`Log? fold change`))
                                      # the case of not choosing top n
                                      if (top_n == 0) {
                                        return (df_filtered)
                                      }

                                      df_top_n <- df_filtered %>%
                                        dplyr::slice(1:top_n) %>%
                                        dplyr::select("Name", `Log? fold change`, `P-value`)
                                      return(df_top_n)

                                    },

                                    # get tail genes by log2 fold change, return df
                                    get_tail_n_genes = function(top_n, log_2_FC, p_value) {
                                      df_filtered <- self$filter_df(log_2_FC, p_value)%>%
                                        arrange(desc(`Log? fold change`))
                                      # the case of not choosing top n
                                      if (top_n == 0) {
                                        return (df_filtered)
                                      }

                                      num_rows <- nrow(df_filtered)
                                      df_tail_n <- df_filtered %>%
                                        dplyr::slice((num_rows - top_n + 1):num_rows) %>%
                                        dplyr::select("Name", `Log? fold change`, `P-value`)
                                      return (df_tail_n)
                                      
                                    },

                                    # get top and tail genes by log2 fold change, return df
                                    get_top_and_tail_n_genes = function(top_n, log_2_FC, p_value){
                                      df_top_n <- self$get_top_n_genes(top_n, log_2_FC, p_value)
                                      df_tail_n <- self$get_tail_n_genes(top_n, log_2_FC, p_value)
                                      
                                      combined_df <- bind_rows(df_top_n, df_tail_n)
                                      return(combined_df)
                                    }
                                    


                                  )
)

DataProcessor <- R6Class("DataProcessor",
                    public = list(
                      diff_expr_file = NULL,
                      gene_folder = NULL,
                      pattern = NULL,

                      
                      initialize = function(diff_expr_file, gene_folder, pattern) {
                        self$diff_expr_file <- diff_expr_file
                        self$gene_folder <- gene_folder
                        self$pattern <- pattern

                      },
                      
                      extract_top_n_gene_feature = function(feature, top_n, pivot, log_2_FC, p_value) {

                        diff_expr <- DifferentialExpressionDataProcessor$new(self$diff_expr_file)
                        print(diff_expr, top_n, log_2_FC)
                        gene_top_n <- diff_expr$get_top_and_tail_n_genes(top_n, log_2_FC, p_value)
                        feature_processor <- GeneFeatureExtractor$new(self$gene_folder, self$pattern, feature)
                        
                        gene_feature_pivot <- feature_processor$extract_feature_pivot(gene_top_n)
                        if (pivot == FALSE) {
                          gene_feature <- feature_processor$extract_feature_unpivot(gene_feature_pivot)
                        } else {
                          gene_feature <- gene_feature_pivot
                        }
                        
                        return(gene_feature)
                      }
                    )
)

