setwd(out_dir)

# preprocess data
cds <- preprocess_cds(cds, num_dim = 50)
cds <- align_cds(cds, alignment_group = "condition")
cds <- reduce_dimension(cds)

# cluster cells
cds <- cluster_cells(cds)
cds <- learn_graph(cds, use_partition = FALSE)

#  reorder cells by pseudotime
cds_order <- order_cells(cds, root_pr_nodes = c("Y_57"))

pdf("UMAP.pdf")
# show all
plot_cells(cds,
           label_cell_groups=TRUE,
           show_trajectory_graph=FALSE,
           color_cells_by = "nowakowski.noglyc",
           group_label_size=5,
           graph_label_size=3) +
  theme(axis.text = element_text(size = 11, face="bold"),
        axis.title = element_text(size = 11, face = "bold"))

mygenotype = levels(factor(cds$condition))
for (g in mygenotype){
  print(plot_cells(cds[,colnames(cds)[cds$condition == g]],
                   label_groups_by_cluster=FALSE,
                   color_cells_by = "nowakowski.noglyc",
                   group_label_size = 5,
                   graph_label_size=3
  ) + annotate("text", x=0, y=10, label= g, size = 7) +
    theme(axis.text = element_text(size = 11, face="bold"),
          axis.title = element_text(size = 11, face = "bold"))
  )
}
dev.off()

pdf("Trajectory.pdf")
plot_cells(cds,
           color_cells_by = "nowakowski.noglyc",
           label_groups_by_cluster=FALSE,
           label_leaves=FALSE,
           label_branch_points=TRUE,
           label_principal_points = FALSE,
           group_label_size=3,
           graph_label_size=3)

for (g in mygenotype){
  print(plot_cells(cds[,colnames(cds)[cds$condition == g]],
                   color_cells_by = "nowakowski.noglyc",
                   label_groups_by_cluster=FALSE,
                   label_leaves=FALSE,
                   label_branch_points=TRUE,
                   label_principal_points = FALSE,
                   group_label_size = 4,
                   graph_label_size=3
  ) + annotate("text", x=10, y=7, label= g, size = 3.5)
  )
}
dev.off()

pdf("Trajectory_pseudotime.pdf")
plot_cells(cds_order,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           group_label_size = 4,
           graph_label_size = 4)

plot_cells(cds_order,
           color_cells_by = "nowakowski.noglyc",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           group_label_size = 4,
           graph_label_size=4)
dev.off()