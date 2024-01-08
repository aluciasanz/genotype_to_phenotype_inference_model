options(stringsAsFactors = F)
library(glmnet)

set.seed(616)

#### variables used 

pos_label = label_allgen[label_allgen > 0]
pos_label_log = log(pos_label)
pos_fof = p_na_b_aa_1st[label_allgen > 0, ]
pos_sof = p_na_b_aa_2nd[label_allgen > 0, ]
pos_cof = cbind(pos_fof, pos_sof)



#### second step linear model on positive eop pairs only with log transform and normality test

#boot_2nd_step_run_mae_res = sapply(1:200, second_step_lm_pos_log, pos_label_log, pos_fof, pos_sof, pos_cof, 1, 'mae')
#boot_2nd_step_run_mse_res = sapply(1:200, second_step_lm_pos_log, pos_label_log, pos_fof, pos_sof, pos_cof, 1, 'mse')

# pdf("p3_f3_lm_training_error_comp.pdf", height = 7.5, width = 6, useDingbats = F)
# training_lm_res = melt(t(boot_2nd_step_run_mae_res[1:4, ]))
# boxplot(value ~ Var2, data = training_lm_res, names = c("null", "FOF", "SOF", "COF"), ylab = "MAE", main = "Training MAE on second step linear model", col = c(rgb(0, 0, 0, 0.5), rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), notch = T)
# dev.off()
# 
# pdf("p3_f3_lm_test_error_comp.pdf", height = 7.5, width = 6, useDingbats = F)
# test_lm_res = melt(t(boot_2nd_step_run_mae_res[5:8, ]))
# boxplot(value ~ Var2, data = test_lm_res, names = c("null", "FOF", "SOF", "COF"), ylab = "MAE", main = "Test set MAE on second step linear model", col = c(rgb(0, 0, 0, 0.5), rgb(1, 0, 0, 0.5), rgb(0, 1, 0, 0.5), rgb(0, 0, 1, 0.5)), notch = T)
# dev.off()


#### non-filtered phage variables

pos_fof_bo_all = p_aa_b_aa_1st[label_allgen > 0, 178:196]
pos_fof_po_all = p_aa_b_aa_1st[label_allgen > 0, 1:177]
pos_fof_all = p_aa_b_aa_1st[label_allgen > 0, ]
pos_sof_all = p_aa_b_aa_2nd[label_allgen > 0, ]
pos_cof_all = cbind(pos_fof_all, pos_sof_all)

#boot_2nd_step_run_mae_res_all = pbsapply(1:200, second_step_lm_pos_log, pos_label_log, pos_fof_all, pos_sof_all, pos_cof_all, 1, 'mae')
#boot_2nd_step_run_mse_res_all = pbsapply(1:200, second_step_lm_pos_log, pos_label_log, pos_fof_all, pos_sof_all, pos_cof_all, 1, 'mse')

step2_meta_results = pblapply(1:200, second_step_lm_pos_log_5model, pos_label_log, pos_fof_bo_all, pos_fof_po_all, pos_fof_all, pos_sof_all, pos_cof_all, 1, 'mae')

boot_2nd_step_run_mae_res_all = sapply(step2_meta_results, function(x) x$stats)

boot_2nd_step_run_mae_res_all_coef_fof_m = sapply(step2_meta_results, function(x) x$coef_trained$coef_fof[, 1])

boot_2nd_step_run_mae_res_all_test_samples_m = rowSums(sapply(step2_meta_results, function(x) x$test_idx_bined))

# training_lm_res_allf = melt(t(boot_2nd_step_run_mae_res_all[1:6, ]))
# 
# par(mar = c(8, 4, 2, 2))
# boxplot(value ~ Var2, data = training_lm_res_allf, names = c("null", "hFOF", "pFOF", "FOF", "SOF", "COF"), ylab = "MAE", main = "Training MAE on second step linear model", col = rgb(0, 0, 0, 0.5), notch = T)
# 
# test_lm_res_allf = melt(t(boot_2nd_step_run_mae_res_all[7:12, ]))
# 
# par(mar = c(8, 4, 2, 2))
# boxplot(value ~ Var2, data = test_lm_res_allf, names = c("null", "hFOF", "pFOF", "FOF", "SOF", "COF"), ylab = "MAE", main = "Test set MAE on second step linear model", col = rgb(0, 0, 0, 0.5), notch = T)
# 

save.image("second_step_all_feature_run_res3_FINAL.RData")
