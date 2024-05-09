library(tidyr)
library(dplyr)
library(readr)
library(readxl)
library(tidyverse)
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

## Get TPM from gene files into a pivot
process_files_tpm_pivot <- function(directory, df, pattern) {
  file_paths <- list.files(directory, pattern = pattern, full.names = TRUE)
  updated_df <- df
  
  # Loop over file paths
  for (file_path in file_paths) {
    # Read the CSV file
    data_df <- read_csv(file_path)
    
    
    # Prepare the filename for the new column
    file_name <- gsub(pattern = " \\(GE\\)\\.csv$", replacement = "", x = basename(file_path))
    file_name_column <- paste(file_name, "TPM", sep = "_")
    
    # Join the data with df based on 'Name' and 'Molecules'
    merged_df <- updated_df %>%
      left_join(data_df, by = c("Name" = "Name")) %>%
      dplyr::select(Name, TPM) %>%
      dplyr::rename(!!file_name_column := TPM)
    
    updated_df[[file_name_column]] <- merged_df[[file_name_column]]
    
    # remove NA
    updated_df <- updated_df %>%
      filter(!if_any(ends_with("TPM"), is.na))
  }
  
  return(updated_df)
}

## unpivot df
unpivot <- function(df_pivot){
  df_unpivot <- df_pivot %>%
    pivot_longer(
      cols = ends_with("_TPM"),
      names_to = "Sample", 
      values_to = "TPM"
    ) %>%
    dplyr::mutate(
      Number = str_extract(Sample, "(?<=_)\\d+(?=_?)"), # Extracts the number with the underscore
      # Number = str_remove(Number, "_"),
      Sample = str_remove(Sample, "_TPM"),  # Removes the '_TPM' suffix
      Sample = str_remove(Sample, "_\\d+_"), # Removes the leading number with underscore
      Gene = str_extract(Sample, "^[^\\-]+"), # Extracts the gene part before the first hyphen
      Sample = str_remove(Sample, "^[^\\-]+\\-"), # Removes the gene part and the following hyphen
      Gender = str_sub(Sample, 1, 1),  # Assumes the gender is denoted by the first character after gene
      Treatment = str_sub(Sample, 3) # Assumes treatment is the rest after removing the gender and hyphen
    ) %>%
    dplyr::select(-Sample)
  print(df_unpivot)
  return(df_unpivot)
}

## Export CSV file from input xls file 
gen_tpm_data <- function(folder, pathways_list, pattern, input_file, output_name) {
  gene_pathway_df <- extract_pathway_gene(input_file, pathways_list)
  
  gene_pathway_tpm <- process_files_tpm(folder, gene_pathway_df, pattern)
  
  write.csv(gene_pathway_tpm, output_name, row.names = FALSE)
}

# data for standard heatmap #######################################
extract_diff_expr_tpm <- function(folder, diff_expr, gene_folder, pattern, top_n, pivot) {
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
  
  gene_tpm_pivot <- process_files_tpm_pivot(folder, combined_df, pattern)
  if (pivot == F) {
    gene_tpm <- unpivot(gene_tpm_pivot)
  } else {
    gene_tpm <- gene_tpm_pivot
  }
  
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
