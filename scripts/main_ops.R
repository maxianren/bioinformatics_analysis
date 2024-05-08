setwd("~/Documents/work/bio_analysis/bioinformatics_analysis")
source("scripts/format_data.R")
source("scripts/draw.R")



# pipe line for heat map & bar chart pathway ####################################################################
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

# pipe line for heat map ####################################################################
## var =========================================================================
gene_folder = "data/LPS/genes"
data_dir = "data/LPS"
out_dir = "out/LPS"
## compare treatment ====
column_level_1_treatment = c("pWT", "pWT","pWT","pWT","pWT","pKO","pKO","pKO","pKO","pKO","pKO")
column_level_2_treatment = c("PBS","PBS","LPS_5_ug","LPS_5_ug","LPS_5_ug", "PBS","PBS","LPS_5_ug","LPS_5_ug","LPS_5_ug","LPS_5_ug")
column_order_treatment = c("_1_pWT-M-PBS_TPM", 
                           "_6_pWT-F-PBS_TPM", 
                           "_3_pWT-M-5_g_TPM", 
                           "_7_pWT-F-5_g_TPM", 
                           "_8_pWT-F-5_g_TPM",
                           "_11_pKO-M-PBS_TPM", 
                           "_15_pKO-F-PBS_TPM", 
                           "_12_pKO-M-5_g_TPM", 
                           "_13_pKO-M-5_g_TPM", 
                           "_16_pKO-F-5_g_TPM", 
                           "_17_pKO-F-5_g_TPM")
column_title_treatment = "Gene Expression - LPS vs PBS, Sorted by Log2 fold change Top 60"

df_diff_expr = extract_diff_expr_tpm(file.path(data_dir,"diff_expr_5_g vs. PBS.csv"), 
                             folder = gene_folder,
                             pattern = "(.*KO.*5_g.*|.*KO.*PBS.*|.*WT.*5_g.*|.*WT.*PBS.*)",
                             top_n = 30,
                             pivot = T)

draw_heatmap_std(df = df_diff_expr, 
                 column_order =column_order_treatment, 
                 column_level_1 = column_level_1_treatment, 
                 column_level_2 = column_level_2_treatment, 
                 column_title = "Gene Expression - LPS vs PBS, Sorted by Log2 fold change Top 60",
                 flag_row_name = T,
                 output = file.path(out_dir, "heatmap_ko_5_g vs. PBS - top60.png"))

####  5g vs 100g vs pbs -----

### compare gene ----
df_diff_expr = extract_diff_expr_tpm(file.path(data_dir,"diff_expr_pWT vs. pKO.csv"), 
                                 folder = gene_folder,
                                 pattern = "(.*KO.*5_g.*|.*KO.*PBS.*|.*WT.*5_g.*|.*WT.*PBS.*)",
                                 top_n = 30,
                                 sort_by = 'Log? fold change')

draw_heatmap_std(df = df_diff_expr, 
                 column_order =column_order_treatment, 
                 column_level_1 = column_level_1_treatment, 
                 column_level_2 = column_level_2_treatment, 
                 column_title = "Gene Expression - WT vs KO, Sorted by Log2 fold change Top 60",
                 flag_row_name = T,
                 output = file.path(out_dir, "heatmap_ko_diff_expr_pWT vs. pKO - top60.png"))


## go analysis ====
log_2_FC = 1
top_n = 10
gene = extract_diff_expr_log2fc(input_file=file.path(data_dir,"diff_expr_5_g vs. PBS.csv"), 
                              log_2_FC=log_2_FC)
draw_go_analysis(gene, 
                 top_n = top_n,
                 title = "Go Analysis: 5_g vs. PBS",
                 output = file.path(out_dir, "go_enrichment_5_g vs. PBS.png"))

# box plot
df_diff_expr_unpivot = extract_diff_expr_tpm(file.path(data_dir,"diff_expr_5_g vs. PBS.csv"), 
                                     folder = gene_folder,
                                     pattern = "(.*KO.*5_g.*|.*KO.*PBS.*|.*WT.*5_g.*|.*WT.*PBS.*)",
                                     top_n = 5,
                                     pivot = F)

draw_box_plot(df_diff_expr_unpivot,
              title = "Box Plot 5_g vs. PBS",
              output = file.path(out_dir, "box_plot_5_g vs. PBS.png"))
