# params
GOI <- c("RBFOX1","SETD1A","SOX2","FOXG1","MAP2","CTIP2","NEUROD2")

# create a subset cds with only GOI
cds_subset <- cds[rowData(cds)$gene_short_name %in% GOI,]

# identify variable genes
gene_fits <- fit_models(cds_subset, model_formula_str = "~condition") # model_formula_str = "~condition + nowakowski_med", etc. 
fit_coefs <- coefficient_table(gene_fits)

# head(fit_coefs)

# make violin plot for each interesting term
my_terms <- fit_coefs %>% filter(term == "conditionFS")
fit_res = my_terms %>% filter (q_value < 0.05) %>%
  select(gene_short_name, term, q_value, estimate)
fit_res_order = arrange(fit_res, q_value)
write.table(fit_res_order,"GOI_fit_res.txt")

pdf("GOI_violin_plots.pdf")
plot_genes_violin(cds_subset, group_cells_by="condition", ncol=2,normalize = FALSE) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
dev.off()

# choose a distribution to base gene expression models on
evaluate_fits(gene_fits)

time_batch_models <- fit_models(cds_subset,
                                model_formula_str = "~condition + nowakowski_med",
                                expression_family="negbinomial")
time_models <- fit_models(cds_subset,
                          model_formula_str = "~condition",
                          expression_family="negbinomial")
# If p_value < 0.05, the factor should be included
compare_models(time_batch_models, time_models) %>% select(gene_short_name, q_value)
