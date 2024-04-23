source("~/result/scripts/format_data.R")
source("~/result/scripts/draw.R")


# pipe line ####################################################################
## var =========================================================================
gene_folder = "~/result/data/gene_files"
data_dir = "~/result/data/"
out_dir = "~/result/out/"
pathways_list_wt = c("Airway Inflammation in Asthma",
                     "Airway Pathology in Chronic Obstructive Pulmonary Disease (COPD)",
                     "Role of IL-17F in Allergic Inflammatory Airway Diseases",
                     "IL-17A Signaling in Airway Cells",
                     "Hepatic Fibrosis / Hepatic Stellate Cell Activation",
                     "Neuroinflammation Signaling Pathway",
                     "Pulmonary Healing Signaling Pathway",
                     "Pulmonary Fibrosis Idiopathic Signaling Pathway",
                     "Small Cell Lung Cancer Signaling",
                     "Role of NFAT in Regulation of the Immune Response",
                     "Hepatic Fibrosis Signaling Pathway"
                     )

pathways_list_ko = c("Airway Pathology in Chronic Obstructive Pulmonary Disease (COPD)",
                     "Airway Inflammation in Asthma",
                     "Role of IL-17F in Allergic Inflammatory Airway Diseases",
                     "Neuroinflammation Signaling Pathway",
                     "IL-17A Signaling in Airway Cells",
                     "Hepatic Fibrosis / Hepatic Stellate Cell Activation",
                     "Pulmonary Healing Signaling Pathway",
                     "Hepatic Fibrosis Signaling Pathway",
                     "Small Cell Lung Cancer Signaling",
                     "Pulmonary Fibrosis Idiopathic Signaling Pathway",
                     "The citric acid (TCA) cycle and respiratory electron transport chain",
                     "Role of NFAT in Regulation of the Immune Response",
                     "Non-Small Cell Lung Cancer Signaling"
                     )
column_level_1_wt_ko = c("PBS", "PBS", "PBS", "PBS", "PBS","HDM", "HDM",  "HDM", "HDM", "HDM")
column_level_2_wt_ko = c("M", "M", "F", "F", "F","M", "M",  "F", "F", "F")

## WT - HDM vs PBS =============================================================
gen_tpm_data(folder = gene_folder,
             pathways_list = pathways_list_wt,  
             pattern = "WT*",
             input_file =  file.path(data_dir, "pathway_wt.xls"),
             output_name =  file.path(data_dir, "/gene_pathway_tpm_wt_all.csv")
)

draw_bar(input_file = file.path(data_dir, "pathway_wt.xls"),
         pathways_list = pathways_list_wt,
         title = "Bar Chart for Canonical Pathways - WT HDM vs PBS",
         output = file.path(out_dir, "bar_chart_wt.png")
)

draw_heatmap(input_file = file.path(data_dir, "/gene_pathway_tpm_wt_all.csv"), 
             column_order = c("_43_WT-M-PBS_TPM", 
                              "_44_WT-M-PBS_TPM", 
                              "_55_WT-F-PBS_TPM", 
                              "_56_WT-F-PBS_TPM",
                              "_57_WT-F-PBS_TPM", 
                              "_46_WT-M-HDM_TPM", 
                              "_48_WT-M-HDM_TPM", 
                              "_58_WT-F-HDM_TPM", 
                              "_60_WT-F-HDM_TPM",
                              "_61_WT-F-HDM_TPM"), 
             column_level_1 = column_level_1_wt_ko, 
             column_level_2 = column_level_2_wt_ko, 
             column_title = "Gene Expression by WT - HDM vs PBS",
             output = file.path(out_dir, "heatmap_wt.png")
             )
## KO - HDM vs PBS =============================================================
gen_tpm_data(folder = gene_folder,
             pathways_list = pathways_list_ko,  
             pattern = "KO*",
             input_file = file.path(data_dir, "pathway_ko.xls"),
             output_name = file.path(data_dir, "gene_pathway_tpm_ko_all.csv")
)

draw_bar(input_file = file.path(data_dir, "pathway_ko.xls"),
         pathways_list = pathways_list_ko,
         title = "Bar Chart for Canonical Pathways - KO HDM vs PBS",
         output = file.path(out_dir, "bar_chart_ko.png")
)

draw_heatmap(input_file = file.path(data_dir, "gene_pathway_tpm_ko_all.csv"), 
             column_order = c("_49_KO-M-PBS_TPM",
                              "_50_KO-M-PBS_TPM",
                              "_62_KO-F-PBS_TPM",
                              "_63_KO-F-PBS_TPM",
                              "_64_KO-F-PBS_TPM",
                              "_72_KO-M-HDM_TPM",
                              "_51_KO-M-HDM_TPM",
                              "_68_KO-F-HDM_TPM",
                              "_69_KO-F-HDM_TPM",
                              "_70_KO-F-HDM_TPM"), 
             column_level_1 = column_level_1_wt_ko, 
             column_level_2 = column_level_2_wt_ko, 
             column_title = "Gene Expression by KO - HDM vs PBS",
             output = file.path(out_dir, "heatmap_ko.png")
)
