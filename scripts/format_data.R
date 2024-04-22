library(tidyr)
library(dplyr)
library(readr)

extract_pathway_gene <- function(file_path, pathways) {
  # Read the CSV file
  df <- read.csv(file_path, , stringsAsFactors = FALSE, check.names = FALSE)
  
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
      left_join(data_df, by = c("Molecule" = "Name")) %>%
      select(Pathway, Molecule, TPM) %>%
      rename(!!file_name_column := TPM)
    
    updated_gene_pathway_df[[file_name_column]] <- merged_df[[file_name_column]]
    
    # remove NA
    updated_gene_pathway_df <- updated_gene_pathway_df %>%
      filter(!if_any(ends_with("TPM"), is.na))
  }
  
  return(updated_gene_pathway_df)
}

gen_tpm_data <- function(folder, pathways_list, pattern, input_file, output_name) {
  gene_pathway_df <- extract_pathway_gene(input_file, pathways_list)
  
  gene_pathway_tpm <- process_files_tpm(folder, gene_pathway_df, pattern)
  
  write.csv(gene_pathway_tpm, output_name, row.names = FALSE)
}


