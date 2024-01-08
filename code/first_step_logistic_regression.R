options(stringsAsFactors = F)

set.seed(106)

# binarize the phenotype and do random forest to distinguish the most important mutations

# phage_na_table = phage_table[which(phage_table$nonsyn_mutation == 1), ]
# rownames(phage_na_table) = paste("phage_mut_nonsyn_", 1:nrow(phage_na_table), sep = "")

#bin_pheno_allgen = as.factor(ifelse(all_labels_raw > 0, 1, 0))

bin_pheno_allgen = as.factor(ifelse(all_labels_raw > 0, 1, 0))
fof_m_bo = p_aa_b_aa_1st[, 178:196] # bacteria mutations
fof_m_po = p_aa_b_aa_1st[, 1:177] # phage mutations including phage ancestor
fof_m_all = p_aa_b_aa_1st
sof_m_all = p_aa_b_aa_2nd
cof_m_all = data.frame(fof_m_all, sof_m_all)

#all_logistic_res_allf = pbsapply(1:200, first_step_logistic_all, bin_pheno_allgen, fof_m_all, sof_m_all, cof_m_all, 1, "class")
step1_meta_results = pblapply(1:200, first_step_logistic_all_5model, bin_pheno_allgen, fof_m_bo, fof_m_po, fof_m_all, sof_m_all, cof_m_all, 1, "class")

all_logistic_res_allf = sapply(step1_meta_results, function(x) x$stats)

all_logistic_res_allf_coef_fof_m = sapply(step1_meta_results, function(x) x$coef_trained$coef_fof[, 1])

all_logistic_res_allf_test_samples_m = rowSums(sapply(step1_meta_results, function(x) x$test_idx_bined))

# train_res_l = melt(t(all_logistic_res_allf[1:6, ]))
# 
# par(las = 2, mar = c(8, 4, 2, 2))
# boxplot(value ~ Var2, data = train_res_l, names = c("Null", "HoMF", "PoMF", "P&HMF", "PCHMF", "CoMF"), ylab = "Classification error", main = "Training error on first step logistic regression", col = rgb(0, 0, 0, 0.5), notch = T)
# 
# test_res_l = melt(t(all_logistic_res_allf[7:12, ]))
# 
# par(las = 2, mar = c(8, 4, 2, 2))
# boxplot(value ~ Var2, data = test_res_l, names = c("Null", "HoMF", "PoMF", "P&HMF", "PCHMF", "CoMF"), ylab = "Classification error", main = "Test error on first step logistic regression", col = rgb(0, 0, 0, 0.5), notch = T)

save.image("first_step_all_feature_logistic_run_res3_FINAL.RData")
