################################
## volcano plot with ggplot2
library(ggplot2 )

# read the input res file
res <- read.csv("D:/Data/setd1a/5-seurat/out/10_01_2023_b4_deg_out/IN-STR_fs3_cntrl3_mast.csv", row.names=1) # input csv for mast
out_file <- "13_01_2023_b4_IN-STR_volcano_plot.pdf"

## BATCH 2 MAST ANALYSIS
# "D:/Data/setd1a/5-seurat/out/10_01_2023_b2_deg_out/EN-V1_fs3_cntrl3_mast.csv"       DONE
# "D:/Data/setd1a/5-seurat/out/10_01_2023_b2_deg_out/RG-early_fs3_cntrl3_mast.csv"    DONE

## BATCH 4 MAST ANALYSIS
# "D:/Data/setd1a/5-seurat/out/10_01_2023_b4_deg_out/EN-V1_fs3_cntrl3_mast.csv"
# "D:/Data/setd1a/5-seurat/out/10_01_2023_b4_deg_out/RG-early_fs3_cntrl3_mast.csv"
# "D:/Data/setd1a/5-seurat/out/10_01_2023_b4_deg_out/IN-STR_fs3_cntrl3_mast.csv"

head(res) 
volcano.input <- data.frame(gene = row.names(res),
                  pvalue = -log10(res$p_val_adj),
                  lfc = res$avg_log2FC)
volcano.input <- na.omit(volcano.input)
volcano.input <- volcano.input[volcano.input$pvalue != Inf,]

# # no_filter ones make more sense in this case
# volcano.input <- data.frame(gene = row.names(res_nofilter),
#                             pvalue = -log10(res_nofilter$padj), 
#                             lfc = res_nofilter$log2FoldChange)
# volcano.input <- na.omit(volcano.input)

# we can run two condition volcano_plot if needed
#volcano_pair1.input <- data.frame(gene = row.names(resMF_pair2),
#                            pvalue = -log10(resMF_pair2$padj), 
#                            lfc = resMF_pair2$log2FoldChange,
#                            baseMean = resMF_pair2$baseMean)

# volcano_pair1.input <- na.omit(volcano_pair1.input)

# volcano.input = volcano_pair1.input

# change subset colors according to specific Genotype
library(dplyr)  # ggplot2 addon for better labels
library(ggrepel) # ggplot2 addon for better labels

# # get miRNA genes
# miRNA_genes <- volcano.input[grepl('MIR.*', volcano.input$gene,perl=TRUE),]$gene

# color column
volcano.input <- volcano.input %>%
  mutate(color = ifelse(volcano.input$lfc >= 0 & volcano.input$pvalue >= 1.3,
                        yes = "Up",
                        no = ifelse(volcano.input$lfc <= 0 & volcano.input$pvalue >= 1.3,
                                    yes = "Down",
                                                no = "Others")))
head(volcano.input)
  # mutate(color = ifelse(volcano.input$gene %in% del_genes$GeneSymbol, 
  #                       yes = "22q11DS_gene", 
  #                       no = ifelse(volcano.input$gene %in% miRNA_genes, 
  #                                   yes = "miRNA", 
  #                                   no = ifelse(volcano.input$gene == "C19ORF63", 
  #                                               yes = "Mirta22/Emc10", 
  #                                               no = "Others"))))

myplot <- ggplot(volcano.input, aes(x = lfc, y = pvalue)) + 
  #  geom_point(aes(colour = factor(color)), size = 2, alpha = 0.8, na.rm = T) + # Make dots bigger
  theme_bw(base_size = 16) + # change theme
  #  ggtitle(label = "Volcano Plot", subtitle = "Simple black & white") + # Add a title
  xlab(expression(log[2] ("FS" / "CT"))) + # x-axis label
  ylab(expression(-log[10]("adjusted p-value"))) + # y-axis label
  #  geom_vline(xintercept = c(-2,2), colour = "darkgrey") + # Add fc cutoffs
  geom_hline(yintercept = 1.3, colour = "darkgrey", linetype=8) + # Add pvalue cutoffs
  geom_vline(xintercept = 0, colour = "black") + # Add 0 lines
  #  scale_colour_gradient(low = "black", high = "black", guide = FALSE) + # Color black
  scale_x_continuous(limits = c(-2.5, 2.5)) # min/max of lfc
#  annotate(geom = "text", 
#         label = "Untreated", 
#         x = -2, y = 85, 
#         size = 7, colour = "black") + # add Untreated text
#  annotate(geom = "text", 
#           label = "Treated", 
#           x = 2, y = 85, 
#           size = 7, colour = "black") + # add Treated text

# Subset table to only show certain gene labels
labelled <- volcano.input[volcano.input$pvalue>15,]
#labelled <- volcano.input[(grepl("MIR.*", volcano.input$gene,perl=TRUE) | volcano.input$gene == "C19ORF63" | volcano.input$gene %in% del_genes$GeneSymbol) & volcano.input$pvalue >1.3,]

# Add layer of color and text annotation to volcano plot.
pdf(out_file)
myplot + 
  scale_colour_manual(values = c("Others" = "grey",
                                 "Up" = "red", 
                                 "Down" = "blue"
)) + scale_y_continuous(trans = "log1p") + 
  geom_point(data = subset(volcano.input, color == 'Others'),
             aes(x = lfc, y = pvalue, color = 'Others'),na.rm = T) +
  geom_point(data = subset(volcano.input, color == 'Up'),
             aes(x = lfc, y = pvalue, color = 'Up')) +
  geom_point(data = subset(volcano.input, color == 'Down'),
             aes(x = lfc, y = pvalue, color = 'Down')) +
  geom_text_repel(data = labelled, 
                  mapping = aes(label = gene), 
                  size = 5,
                  fontface = 'bold', 
                  color = 'black',
                  box.padding = unit(0.5, "lines"),
                  point.padding = unit(0.5, "lines"))
dev.off()
