# params
cell_type <- "RG-early"

# reload and preprocess raw dataset
cds <- load_monocle_objects(directory_path = 'D:/Data/setd1a/6-monocle3/out/cds_object.rds')
cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds)
cds <- cluster_cells(cds, resolution=1e-5)

# Graph-autocorrelation analysis using predefined cell type
COI_cds <- cds[,grepl("RG-early", colData(cds)$nowakowski_med, ignore.case=TRUE)]

### show the user defined cell population on the trajectory map
# NOTE: This is a very useful function to determine if user defined cell population is homogeneous!!!
plot_cells(COI_cds,
           color_cells_by="cluster",
           label_cell_groups=TRUE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           group_label_size = 6,
           graph_label_size=4
)
### Rank genes with knn. This step will take ~15 minutes!!
pr_graph_test_res <- graph_test(COI_cds, neighbor_graph="knn", cores=12)
pr_graph_test_res_order <-pr_graph_test_res[order(pr_graph_test_res$morans_I,decreasing = TRUE),]
pr_graph_test_res_order

### identify cell type specific driver genes
pr_deg_ids <- row.names(subset(pr_graph_test_res_order, q_value < 10^-40 & morans_I > 0.3))

### group these genes into modules
gene_module_df <- find_gene_modules(COI_cds[pr_deg_ids,], resolution=0.008)

## group cells
cell_group_df <- tibble::tibble(cell=row.names(colData(COI_cds)),
                                cell_group=partitions(cds)[colnames(COI_cds)])
agg_mat <- aggregate_gene_expression(COI_cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))
colnames(agg_mat) <- stringr::str_c("Partition ", colnames(agg_mat))

pheatmap::pheatmap(agg_mat, cluster_rows=TRUE, cluster_cols=TRUE,
                   scale="column", clustering_method="ward.D2",
                   fontsize=6)
## plot top modules
plot_cells(COI_cds,
           genes=gene_module_df %>% filter(module %in% c(8, 36, 31, 29)),
           group_cells_by="partition",
           color_cells_by="partition",
           show_trajectory_graph=FALSE)

#########################################
## Finding genes that change as a function of pseudotime
### find the genes that are differentially expressed on the different paths through the trajectory?
###. Clustering
cds <- cluster_cells(cds)
plot_cells(cds,
           color_cells_by = "partition",
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_principal_points = FALSE,
           group_label_size = 6,
           graph_label_size=3
)
##. Learn the trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

### test whether cells at similar positions on the trajectory have correlated expression
### If yes, what are the driver genes for a specific position
### The morans_I [-1, +1]. A value of 0 indicates no effect, while +1 indicates perfect positive autocorrelation
### and suggests that nearby cells have very similar values of a gene's expression.
ciliated_cds_pr_test_res <- graph_test(cds, neighbor_graph="principal_graph", cores=4)
ciliated_cds_pr_test_res_order <- ciliated_cds_pr_test_res[order(ciliated_cds_pr_test_res$morans_I,decreasing = TRUE),]
pr_deg_ids <- row.names(subset(ciliated_cds_pr_test_res_order, q_value < 0.05 & n_cells > 100 & morans_I > 0.3))

### gene UMAP clustering
gene_module_df <- find_gene_modules(cds[pr_deg_ids,], resolution=c(10^seq(-6,-1)))

### cell clustering
cell_group_df <- tibble::tibble(cell=row.names(colData(cds)),
                                cell_group=partitions(cds)[colnames(cds)])

agg_mat <- aggregate_gene_expression(cds, gene_module_df, cell_group_df)
row.names(agg_mat) <- stringr::str_c("Module ", row.names(agg_mat))

## show differential pattern
pheatmap::pheatmap(agg_mat,
                   scale="column", clustering_method="ward.D2")

## select gene modules according to the heatmap
plot_cells(cds,
           genes=gene_module_df %>% filter(module %in% c(1, 11,4)),
           color_cells_by = "nowakowski_med",
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

# extract driver genes for one or more interested modules
genes=gene_module_df %>% filter(module %in% c(8))
# a gene's dynamics along a single path.

## get GOI symbols
myGOI <-  c("NEUROD2")#head(row.names(ciliated_cds_pr_test_res1), n=6) #c("PAX6", "SOX2", "MAP2","CTIP2","SETD1A")

## define cell lineage
lineage_cds <- cds[rowData(cds)$gene_short_name %in% myGOI,]
lineage_cds_order <- order_cells(lineage_cds, root_pr_nodes = c("Y_57"))


## plot genes dynamics as a function of pseudotime
plot_genes_in_pseudotime(lineage_cds_order,
                         color_cells_by="pseudotime",
                         min_expr=0.01)

p1 <- plot_genes_in_pseudotime(lineage_cds_order[,colnames(lineage_cds_order)[lineage_cds_order$condition == 'CT']],
                               color_cells_by="pseudotime",
                               min_expr=0.05) + annotate("text", x=0, y=30000, label= 'CT', size = 3.5)

p2 <- plot_genes_in_pseudotime(lineage_cds_order[,colnames(lineage_cds_order)[lineage_cds_order$condition == 'FS']],
                               color_cells_by="pseudotime",
                               min_expr=0.05) + annotate("text", x=0, y=30000, label= 'FS', size = 3.5)
p1 + p2

## plot genes on trajectory map
plot_cells(cds,
           genes=head(genes$id, n=10),
           color_cells_by = "nowakowski_med",
           label_cell_groups=FALSE,
           show_trajectory_graph=TRUE,
)
lineage_cds <- cds[rowData(cds)$gene_short_name %in% c("MAP2"),]
for (g in mygenotype){
  lineage_cds_geno = lineage_cds[,colnames(lineage_cds)[lineage_cds$condition == g]]
  lineage_cds_geno_order <- order_cells(lineage_cds_geno, root_pr_nodes = c("Y_35"))
  print(plot_genes_in_pseudotime(lineage_cds_geno_order,
                                 min_expr=0.05
  ) + annotate("text", x=0, y=30000, label= g, size = 3.5)
  )
}
print(plot_cells(lineage_cds_geno_order,
                 color_cells_by = "nowakowski_med",
                 label_groups_by_cluster=FALSE,
                 label_leaves=FALSE,
                 label_branch_points=FALSE,
                 label_principal_points = FALSE,
                 group_label_size = 4,
                 graph_label_size=3
)
)

# Analyzing subset of cell population in single-cell trajectories
### select the banches
cds_subset <- choose_cells(cds)

###  identify genes with interesting patterns of expression that fall only within the region of the trajectory you selected
subset_pr_test_res <- graph_test(cds_subset, neighbor_graph="principal_graph", cores=4)
pr_deg_ids <- row.names(subset(subset_pr_test_res, q_value < 0.05))

### Grouping these genes into modules can reveal fate specific genes or those that are activate immediate prior to or following the branch point:
gene_module_df <- find_gene_modules(cds_subset[pr_deg_ids,], resolution=0.001)

### organize the modules by their similarity (using hclust) over the trajectory
agg_mat <- aggregate_gene_expression(cds_subset, gene_module_df)
module_dendro <- hclust(dist(agg_mat))
gene_module_df$module <- factor(gene_module_df$module,
                                levels = row.names(agg_mat)[module_dendro$order])

plot_cells(cds_subset,
           genes=gene_module_df,
           label_cell_groups=FALSE,
           show_trajectory_graph=FALSE)

## Save monocle objects
### The Monocle3 cell_data_set includes the UMAP models and nearest neighbor indexes, which are not R objects
#save_monocle_objects(cds=cds, directory_path='E:/recovered/Bin_folder/shp_work/CTN/projects/2_setd1a/scRNAseq/monocle3_final', comment='This is setd1a CT3 FS3 scanpy processed. Stored 2022-10-24.')
#cds <- load_monocle_objects(directory_path='E:/recovered/Bin_folder/shp_work/CTN/projects/2_setd1a/scRNAseq/monocle3')


