library(tidyr)
library(dplyr)
library(readr)
library(readxl)
# data for pathway-gene heatmap ####
## Extract genes for pathway
extract_pathway_gene <- function(file_path, pathways) {
  # Read the xls file
  df <- readxl::read_excel(file_path, skip = 1, col_names = TRUE)
  
  # Split the 'Molecules' column into multiple rows
  long_df <- df %>%
    separate_rows(Molecules, sep = ",") %>%
    filter(`Ingenuity Canonical Pathways` %in% pathways)
  
  # Select and rename the columns
  final_df <- long_df %>%
    select(Pathway = `Ingenuity Canonical Pathways`, Molecule = Molecules)

  final_df$Molecule <- tolower(final_df$Molecule)
  final_df$Molecule <- sub("^(.)", "\\U\\1", final_df$Molecule, perl = TRUE)

  return(final_df)
}

## Get TPM for each pathway's genes
process_files_tpm <- function(directory, gene_pathway_df, pattern) {
  file_paths <- list.files(directory, pattern = pattern, full.names = TRUE)
  updated_gene_pathway_df <- gene_pathway_df
  
  # Loop over file paths
  for (file_path in file_paths) {
    # Read the CSV file
    data_df <- read_csv(file_path)
    
    
    # Prepare the filename for the new column
    file_name <- gsub(pattern = " \\(GE\\)\\.csv$", replacement = "", x = basename(file_path))
    file_name_column <- paste(file_name, "TPM", sep = "_")
    
    # Join the data with gene_pathway_df based on 'Name' and 'Molecules'
    merged_df <- updated_gene_pathway_df %>%
      left_join(data_df, by = c("Name" = "Name")) %>%
      dplyr::select(Name, TPM) %>%
      dplyr::rename(!!file_name_column := TPM)
    
    updated_gene_pathway_df[[file_name_column]] <- merged_df[[file_name_column]]
    
    # remove NA
    updated_gene_pathway_df <- updated_gene_pathway_df %>%
      filter(!if_any(ends_with("TPM"), is.na))
  }
  
  return(updated_gene_pathway_df)
}
## Export CSV file from input xls file 
gen_tpm_data <- function(folder, pathways_list, pattern, input_file, output_name) {
  gene_pathway_df <- extract_pathway_gene(input_file, pathways_list)
  
  gene_pathway_tpm <- process_files_tpm(folder, gene_pathway_df, pattern)
  
  write.csv(gene_pathway_tpm, output_name, row.names = FALSE)
}

# data for standard heatmap #######################################
extract_diff_expr_tpm <- function(folder, diff_expr, gene_folder, pattern, top_n, sort_by) {
  df_diff_expr <- read.csv(diff_expr, check.names = FALSE)
  
  df_sort <- df_diff_expr %>%
    filter(!is.nan(`P-value`)) %>%
    arrange(desc(`Log? fold change`))
  
  # select top N gene by 
  df_top_n <- df_sort %>%
    dplyr::slice(1:top_n) %>%
    dplyr::select("Name", `Log? fold change`, `P-value`)
  # select tail N gene by 
  num_rows <- nrow(df_sort)
  df_tail_n <- df_sort %>%
    dplyr::slice((num_rows - top_n + 1):num_rows) %>%
    dplyr::select("Name", `Log? fold change`, `P-value`)
  
  combined_df <- bind_rows(df_top_n, df_tail_n)
  
  print(combined_df)

 

  gene_tpm <- process_files_tpm(folder, combined_df, pattern)
  
  return(gene_tpm)
  
}

extract_diff_expr_log2fc <- function(input_file, log_2_FC){
  df_diff_expr <- read.csv(input_file, check.names = FALSE)
  
  # filter according to log_2_FC
  df <- df_diff_expr %>%
    filter(!is.na(`P-value`) & `Log? fold change` > log_2_FC) %>%
    arrange(desc(`Log? fold change`))
    
  df_top <- df %>%
    # dplyr::slice(1:top_n) %>%
    dplyr::select('Name')

  # get name of genes as vector
  name_vector <- df_top %>%
    pull(Name) 
  return(name_vector)
}
