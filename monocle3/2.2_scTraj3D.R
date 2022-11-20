# using 3D trajectories
cds_sub <- preprocess_cds(cds_sub, num_dim = 50)
cds_3d <- reduce_dimension(cds_sub, max_components = 3)
cds_3d <- cluster_cells(cds_3d)
cds_3d <- learn_graph(cds_3d)
cds_3d <- order_cells(cds_3d, root_pr_nodes = c("Y_35"))

plot_cells_3d(cds_3d,
              color_cells_by = "nowakowski_med",
)