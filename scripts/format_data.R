library(R6)
library(tidyr)
library(dplyr)
library(readr)
library(readxl)
library(tidyverse)

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
                           }
                         )
)

TPMExtractor <- R6Class("TPMExtractor",
                        inherit = GeneDataProcessor,
                        
                        public = list(
                          extract_tpm_pivot = function(df) {
                            file_paths <- self$get_file_paths()
                            updated_df <- df
                            
                            for (file_path in file_paths) {
                              data_df <- read_csv(file_path)
                              file_name <- gsub(pattern = " \\(GE\\)\\.csv$", replacement = "", x = basename(file_path))
                              file_name_column <- paste(file_name, "TPM", sep = "_")
                              
                              merged_df <- updated_df %>%
                                left_join(data_df, by = c("Name" = "Name")) %>%
                                dplyr::select(Name, TPM) %>%
                                dplyr::rename(!!file_name_column := TPM)
                              
                              updated_df[[file_name_column]] <- merged_df[[file_name_column]]
                              updated_df <- updated_df %>%
                                filter(!if_any(ends_with("TPM"), is.na))
                            }
                            return(updated_df)
                          },
                          
                          extract_tpm_unpivot = function(df_pivot) {
                            df_unpivot <- df_pivot %>%
                              pivot_longer(
                                cols = ends_with("_TPM"),
                                names_to = "Sample", 
                                values_to = "TPM"
                              ) %>%
                              dplyr::mutate(
                                Number = str_extract(Sample, "(?<=_)\\d+(?=_?)"),
                                Sample = str_remove(Sample, "_TPM"),
                                Sample = str_remove(Sample, "_\\d+_"),
                                Gene = str_extract(Sample, "^[^\\-]+"),
                                Sample = str_remove(Sample, "^[^\\-]+\\-"),
                                Gender = str_sub(Sample, 1, 1),
                                Treatment = str_sub(Sample, 3)
                              ) %>%
                              dplyr::select(-Sample)
                            return(df_unpivot)
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
                                    
                                    get_top_n_genes = function(top_n) {
                                      df_sort <- self$df_diff_expr %>%
                                        filter(!is.nan(`P-value`)) %>%
                                        arrange(desc(`Log? fold change`))
                                      
                                      df_top_n <- df_sort %>%
                                        dplyr::slice(1:top_n) %>%
                                        dplyr::select("Name", `Log? fold change`, `P-value`)
                                      
                                      num_rows <- nrow(df_sort)
                                      df_tail_n <- df_sort %>%
                                        dplyr::slice((num_rows - top_n + 1):num_rows) %>%
                                        dplyr::select("Name", `Log? fold change`, `P-value`)
                                      
                                      combined_df <- bind_rows(df_top_n, df_tail_n)
                                      return(combined_df)
                                    },
                                    
                                    filter_by_log2_fc = function(log_2_FC) {
                                      df_filtered <- self$df_diff_expr %>%
                                        filter(!is.na(`P-value`) & `Log? fold change` > log_2_FC) %>%
                                        arrange(desc(`Log? fold change`)) %>%
                                        dplyr::select('Name')
                                      
                                      name_vector <- df_filtered %>%
                                        pull(Name)
                                      return(name_vector)
                                    }
                                  )
)

DataProcessor <- R6Class("DataProcessor",
                    public = list(
                      diff_expr_file = NULL,
                      gene_folder = NULL,
                      pattern = NULL,
                      top_n = NULL,
                      pivot = NULL,
                      
                      initialize = function(diff_expr_file, gene_folder, pattern, top_n, pivot = TRUE) {
                        self$diff_expr_file <- diff_expr_file
                        self$gene_folder <- gene_folder
                        self$pattern <- pattern
                        self$top_n <- top_n
                        self$pivot <- pivot
                      },
                      
                      extract_diff_expr_tpm = function() {
                        diff_expr <- DifferentialExpressionDataProcessor$new(self$diff_expr_file)
                        combined_df <- diff_expr$get_top_n_genes(self$top_n)
                        tpm_processor <- TPMExtractor$new(self$gene_folder, self$pattern)
                        
                        gene_tpm_pivot <- tpm_processor$extract_tpm_pivot(combined_df)
                        if (self$pivot == FALSE) {
                          gene_tpm <- tpm_processor$extract_tpm_unpivot(gene_tpm_pivot)
                        } else {
                          gene_tpm <- gene_tpm_pivot
                        }
                        
                        return(gene_tpm)
                      }
                    )
)
