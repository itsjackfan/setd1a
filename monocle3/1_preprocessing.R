library(monocle3)

# params
exp_in <- 'D:/Data/setd1a/6-monocle3/in/data_mat_final_mon.csv'
cm_in <- 'D:/Data/setd1a/6-monocle3/in/obs_mon.csv'
gm_in <- 'D:/Data/setd1a/6-monocle3/in/var_mon.csv'

out_dir <- 'D:/Data/setd1a/6-monocle3/out'

# read in data
expression_matrix <- read.table(
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  exp_in,
  row.names = 1
)

cell_metadata <- read.table(
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  cm_in,
  row.names = 1
)

gene_metadata <- read.table(
  sep = ',',
  header = TRUE,
  check.names = FALSE,
  gm_in,
  row.names = 1
)
gene_metadata$gene_short_name = row.names(gene_metadata)

# create CDS object
cds <- new_cell_data_set(t(expression_matrix),
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_metadata)

# Optionally, you can save the object and load it back.
save_monocle_objects(cds = cds, 
                     directory_path = out_dir, comment='This is setd1a CT3-FS3 raw. Stored 14-11-2022.')
cds <- load_monocle_objects(directory_path = out_dir)

# generate reference Seurat/sData object
sdata <- CreateSeuratObject(
  t(expression_matrix),
  project = 'FULL',
  assay = 'RNA',
  meta.data = cell_metadata
)
