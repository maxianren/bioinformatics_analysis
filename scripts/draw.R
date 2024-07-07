source("scripts/format_data.R")
library(stringr)
library(R6)
library(ggplot2)
library(tidyr)
library(dplyr)
library(ComplexHeatmap)
library(readxl)
library(clusterProfiler)
library(org.Mm.eg.db)
library(AnnotationDbi)
library(enrichplot)
library(gridExtra)
library(corrplot)
library(circlize)

PlottingBase <- R6Class("PlottingBase",
                        public = list(
                          data_dir = NULL,
                          gene_folder = NULL,
                          output = NULL,
                          width = NULL,
                          height = NULL,
                          res = NULL,
                          bg = NULL,

                          initialize = function(data_dir, gene_folder, output, width, height, res, bg) {
                            self$data_dir <- data_dir
                            self$gene_folder <- gene_folder
                            self$output <- output
                            self$width <- width
                            self$height <- height
                            self$res <- res
                            self$bg <- bg
                          },

                          save_plot = function(plot) {
                            png(self$output, width = self$width, height = self$height, res = self$res, bg = self$bg)
                            print(plot)
                            dev.off()
                          }
                        )
)

HeatmapPlot <- R6Class("HeatmapPlot",
                       inherit = PlottingBase,

                       public = list(
                         draw = function(diff_expr_file, column_title, flag_row_name, top_n, pattern) {
                           # var
                           feature = "TPM"

                           df_gene <- GeneFeatureExtractor$new(self$gene_folder, pattern, feature)$get_feature_hierarchy()

                           column_order <- df_gene$file_name_clean
                           column_level_1 <- df_gene$treatment
                           column_level_2 <- df_gene$genotype
                           column_level_3 <- df_gene$gender

                           column_annotation <- HeatmapAnnotation(
                             Treatment = column_level_1,
                             Genotype = column_level_2,
                             Gender = column_level_3,
                             annotation_name_side = "right"
                           )
                           # get dataframe
                           data <- DataProcessor$new(diff_expr_file = file.path(self$data_dir,diff_expr_file),
                                                     gene_folder = gene_folder,
                                                     pattern = pattern)
                           df_diff_expr <- data$extract_top_and_tail_n_gene_feature(top_n = top_n,
                                                                           pivot = T,
                                                                           feature = feature,
                                                                           log_2_FC = 2, 
                                                                           p_value = 0.05)
                           # plot
                           data_matrix <- df_diff_expr %>%
                             dplyr::select(all_of(column_order)) %>%
                             as.matrix()

                           rownames(data_matrix) <- df_diff_expr$Name
                           
                           # max_limit = 4
                           # min_limit = -4
                           
                           z_score_matrix <- t(scale(t(data_matrix)))
                           # z_score_matrix <- pmin(pmax(z_score_matrix, min_limit), max_limit)


                           p <- Heatmap(
                             z_score_matrix,
                             name = "Expression",
                             col = colorRampPalette(c("blue", "white", "red"))(75),
                             cluster_rows = TRUE,
                             cluster_columns = FALSE,
                             show_row_names = flag_row_name,
                             top_annotation = column_annotation,
                             column_title = column_title,
                             row_title = "Gene",
                             row_title_side = "left",
                             row_names_gp = gpar(fontsize = 18),
                             column_names_gp = gpar(fontsize = 18),
                             heatmap_legend_param = list(title = "Expression Level"),
                             column_names_rot = 45,
                             row_names_side = "left",
                             column_names_side = "bottom"
                           )
                           # save the image
                           self$save_plot(p)
                         }
                       )
)


GoAnalysisPlot <- R6Class("GoAnalysisPlot",
                          inherit = PlottingBase,
                          
                          public = list(
                            
                            drawDotPlot = function(diff_expr_file, top_n, title) {
                              # var
                              log_2_FC = 2
                              p_value = 0.05
                              # get dataframe
                              data = DifferentialExpressionDataProcessor$new(file_path = file.path(data_dir,diff_expr_file))
                              genes = data$filter_gene(log_2_FC, p_value)
                              
                              # plot
                              GO_results <- enrichGO(gene = genes, OrgDb = "org.Mm.eg.db", keyType = "SYMBOL", ont = 'ALL')
                              plot <- dotplot(GO_results, showCategory = top_n, split = "ONTOLOGY", title = title) + facet_grid(ONTOLOGY ~ ., scale = "free")

                              # save the image
                              self$save_plot(plot)
                            }
                          )
)

KEGGPlot <- R6Class("KEGGAnalysisPlot",
                          inherit = PlottingBase,
                          
                          private = list(
                            
                            replace_ids_with_symbols = function(gene_ids, mapping_df) {
                              # Split the gene IDs on '/'
                              split_ids <- strsplit(gene_ids, "/")
                              
                              # Replace each ID with its corresponding symbol
                              symbols <- sapply(split_ids, function(ids) {
                                symbols <- mapping_df$SYMBOL[match(ids, mapping_df$ENTREZID)]
                                paste(symbols, collapse="/")
                              })
                              
                              return(symbols)
                            }
                          ),
                          
                          public = list(

                            # KEGG Analysis
                            keggAnalysis = function(diff_expr_file, top_n) {

                              # get dataframe
                              data = DifferentialExpressionDataProcessor$new(file_path = file.path(data_dir,diff_expr_file))
                              genes = data$get_filtered_genes_entrez_id(log_2_FC = 2, p_value = 0.05, org.Mm.eg.db)
                              # genes = data$filter_gene(log_2_FC, p_value)

                              KEGG_results <- enrichKEGG(
                                gene=genes$ENTREZID,
                                organism = "mmu",
                                keyType = "kegg",
                                pvalueCutoff = 0.05,
                                pAdjustMethod = "BH",
                              )
                              return(KEGG_results)

                            },
                            #KEGG dot plot
                            drawDotPlot = function(diff_expr_file, top_n, title){
                              
                              KEGG_results <- self$keggAnalysis(diff_expr_file, top_n)
                              plot <- dotplot(KEGG_results, showCategory = top_n, title = title)
                              # save the image
                              self$save_plot(plot)
                            },
                            # KEGG cnet plot
                            drawCNetPlot = function(diff_expr_file, top_n, title){
                              data = DifferentialExpressionDataProcessor$new(file_path = file.path(data_dir,diff_expr_file))
                              genes = data$get_filtered_genes_entrez_id(log_2_FC = 2, p_value = 0.05, org.Mm.eg.db)

                              KEGG_results = self$keggAnalysis(diff_expr_file, top_n)
                              
                              #replace ids with gene name for better readability
                              KEGG_results@result$geneID <- private$replace_ids_with_symbols(KEGG_results@result$geneID, genes)


                              plot <- cnetplot(KEGG_results,
                                               showCategory = top_n,
                                               categorySize = "pvalue",
                                               # foldChange=genes,
                                               circular = TRUE,
                                               color.params = list(edge = TRUE),
                                               node_label = 'all',
                                               edge_label = FALSE) + 
                                ggtitle(title) +
                                theme(plot.title = element_text(hjust = 0.5, size = 20),
                                      text = element_text(size = 14))
                              # save the image
                              self$save_plot(plot)
                            }
                          )

)

BoxPlot <- R6Class("BoxPlot",
                   inherit = PlottingBase,
                   
                   public = list(
                     draw = function(diff_expr_file, title, pattern) {
                       # var
                       feature = "TPM"

                       # process data
                       data <- DataProcessor$new(diff_expr_file = file.path(self$data_dir,diff_expr_file), 
                                                 gene_folder = self$gene_folder, 
                                                 pattern = pattern)
                       df_diff_expr_unpivot <- data$extract_top_and_tail_n_gene_feature(top_n = top_n, 
                                                                               pivot = F,
                                                                               feature = feature, 
                                                                               log_2_FC = 2, 
                                                                               p_value = 0.05)
                       # plot
                       plot <- ggplot(df_diff_expr_unpivot, aes(x = Name, y = TPM, fill = Treatment)) +
                         geom_boxplot() +
                         labs(title = title, x = "Gene", y = "TPM") +
                         theme_minimal()

                       # save the image
                       self$save_plot(plot)
                     }
                   )
)

CorrelationPlot <- R6Class("BoxPlot",
                   inherit = PlottingBase,
                   
                   public = list(
                     # correlation Analysis
                     correlationAnalysis = function(diff_expr_file, top_n) {
                       
                       data <- DataProcessor$new(diff_expr_file = file.path(self$data_dir,diff_expr_file), 
                                                 gene_folder = self$gene_folder, 
                                                 pattern = pattern)
                       df_diff_expr_unpivot <- data$extract_top_n_gene_feature(top_n = top_n, 
                                                                               pivot = T,
                                                                               feature = feature, 
                                                                               log_2_FC = 2, 
                                                                               p_value = 0.05)
                       # correlation matrix
                       df_matrix <- as.matrix(df_diff_expr_unpivot[, -c(1, 2, 3)])
                       rownames(df_matrix) <- df_diff_expr_unpivot$Name

                       aggregated_correlation <- cor(t(df_matrix), method="pearson")
                       
                       return(aggregated_correlation)
                       
                     },
                     drawMatrixPlot = function(diff_expr_file, title, pattern) {
                       # correlation analysis
                       corr_results <- self$correlationAnalysis(diff_expr_file, top_n)

                       # plot
                       sector.names <- colnames(corr_results)
                       
                       # self$save_plot(plot) doesn't work
                       png(filename = self$output, 
                           width = self$width, 
                           height = self$height, 
                           res = self$res, 
                           bg = self$bg)

                       corrplot(corr_results, method = "color", type = "upper", order = "hclust",
                                addCoef.col = "black",
                                tl.col = "black",
                                tl.srt = 45,
                                cl.pos = "b",
                                cl.cex = 0.75)
                       title(main = title, line = +2)

                       dev.off()

                       # # save the image
                       # self$save_plot(plot)
                     },
                     drawChordPlot = function(diff_expr_file, title, pattern) {
                       # correlation analysis
                       corr_results <- self$correlationAnalysis(diff_expr_file, top_n)

                       # self$save_plot(plot) doesn't work
                       png(filename = self$output, 
                           width = self$width, 
                           height = self$height, 
                           res = self$res, 
                           bg = self$bg)

                       chordDiagram(corr_results, 
                                    transparency = 0.3)
                       title(main = title, line = 0)

                       dev.off()
                       
                       # # save the image
                       # self$save_plot(plot)
                     }
                   )
)