options(stringsAsFactors = F)
library(ridge)
library(Matrix)
library(gridExtra)
library(glmnet)
set.seed(616)
#hmcol2 = c("red", "grey", "grey", "green", "green", "green", "green", "green")

##dist on eop sup figure

# sptest = shapiro.test(log(label_allgen[label_allgen > 0]))
# log_pos_label = log(label_allgen[label_allgen > 0])
# 
# pdf("p3_sf1_eop_distribution.pdf", height = 10, width = 10, useDingbats = F)
# par(mfrow = c(2,2))
# hist(label_allgen, breaks = 100, col = 'gray', main = "Distribution of all EOP values", xlab = "EOP value")
# hist(label_allgen[label_allgen > 0], breaks = 40, col = 'gray', main = "Distribution of positive EOP values only", xlab = "EOP value")
# hist(log_pos_label, breaks = 40, col = 'gray', main = "Distribution of log positive EOP values", xlab = "log(EOP)")
# dev.off()
# 
# pdf("p3_sf1_qqplot.pdf", height = 5, width = 5, useDingbats = F)
# ggqqplot(log_pos_label, main = "Q-Q plot for the log positive EOP values")
# dev.off()

# load("first_step_all_feature_logistic_run_res3.RData")
# load("second_step_all_feature_run_res3.RData")

# Get the median and average of the logistic results for each matrix
#1st lambda from bootstrap
mean_1st_fof_bo_lambda = mean(all_logistic_res_allf[13, ])#fof_bo_lambda
mean_1st_fof_po_lambda = mean(all_logistic_res_allf[14, ])# "fof_po_lambda"
mean_1st_fof_lambda = mean(all_logistic_res_allf[15, ])#fof_lamba
mean_1st_sof_lambda = mean(all_logistic_res_allf[16, ])#sof_lambda
mean_1st_cof_lambda = mean(all_logistic_res_allf[17, ])#cof_lambda
# 
median_1st_fof_bo_lambda = median(all_logistic_res_allf[13, ])
median_1st_fof_po_lambda = median(all_logistic_res_allf[14, ])
median_1st_fof_lambda = median(all_logistic_res_allf[15, ])
median_1st_sof_lambda = median(all_logistic_res_allf[16, ])
median_1st_cof_lambda = median(all_logistic_res_allf[17, ])

#2nd lambda from bootstrap
mean_2nd_fof_bo_lambda = mean(boot_2nd_step_run_mae_res_all[13, ])
mean_2nd_fof_po_lambda = mean(boot_2nd_step_run_mae_res_all[14, ])
mean_2nd_fof_lambda = mean(boot_2nd_step_run_mae_res_all[15, ])
mean_2nd_sof_lambda = mean(boot_2nd_step_run_mae_res_all[16, ])
mean_2nd_cof_lambda = mean(boot_2nd_step_run_mae_res_all[17, ])

median_2nd_fof_bo_lambda = median(boot_2nd_step_run_mae_res_all[13, ])
median_2nd_fof_po_lambda = median(boot_2nd_step_run_mae_res_all[14, ])
median_2nd_fof_lambda = median(boot_2nd_step_run_mae_res_all[15, ])
median_2nd_sof_lambda = median(boot_2nd_step_run_mae_res_all[16, ])
median_2nd_cof_lambda = median(boot_2nd_step_run_mae_res_all[17, ])

####################################
########  final model  #############
####################################
#data
phage_mut_cnt = dim(phage_only_eqtl_m_allmut_allgen)[1]
bac_mut_cnt = dim(bac_only_eqtl_m_allmut_allgen)[1]
phage_cnt = dim(phage_only_eqtl_m_allmut_allgen)[2]
bac_cnt = dim(bac_only_eqtl_m_allmut_allgen)[2]
phage_label = c("phage_an", paste0("phage_mut_", 1:(phage_mut_cnt-1))) #not used
bac_label = c("bac_an", paste0("bac_mut_", 1:(bac_mut_cnt-1))) #not used

all_label = as.factor(ifelse(label_allgen > 0, 1, 0))
all_fof_bo = p_aa_b_aa_1st[, 178:196]
all_fof_po = p_aa_b_aa_1st[, 1:177]
all_fof = p_aa_b_aa_1st
all_sof = p_aa_b_aa_2nd
all_cof = cbind(all_fof, all_sof)

#step1 logistic regression with fof bac only
final_lr_on_all_fof_bo = glmnet(all_fof_bo, all_label, family = "binomial", alpha = 1, lambda = mean_1st_fof_bo_lambda, standardize = F)

first_step_fof_bo_yhat = predict(final_lr_on_all_fof_bo, newx = all_fof_bo, type = "class")
first_step_fof_bo_alpha = final_lr_on_all_fof_bo$a0
first_step_fof_bo_bac_beta = as.numeric(final_lr_on_all_fof_bo$beta)

first_step_fof_bo_bac_beta_l = data.frame(x = 1:bac_mut_cnt, y = 1, coef = first_step_fof_bo_bac_beta)

final_1st_fof_bo_min = min(final_lr_on_all_fof_bo$beta)
final_1st_fof_bo_max = max(final_lr_on_all_fof_bo$beta)

#step1 logistic regression with fof phage only
final_lr_on_all_fof_po = glmnet(all_fof_po, all_label, family = "binomial", alpha = 1, lambda = mean_1st_fof_po_lambda, standardize = F)

first_step_fof_po_yhat = predict(final_lr_on_all_fof_po, newx = all_fof_po, type = "class")
first_step_fof_po_alpha = final_lr_on_all_fof_po$a0
first_step_fof_po_phage_beta = as.numeric(final_lr_on_all_fof_po$beta)

first_step_fof_po_phage_beta_l = data.frame(x = 1:phage_mut_cnt, y = 1, coef = first_step_fof_po_phage_beta)

final_1st_fof_po_min = min(final_lr_on_all_fof_po$beta)
final_1st_fof_po_max = max(final_lr_on_all_fof_po$beta)

#step1 logistic regression with fof (linear model with phage and bacteria)
final_lr_on_all_fof = glmnet(all_fof, all_label, family = "binomial", alpha = 1, lambda = mean_1st_fof_lambda, standardize = F)

first_step_fof_yhat = predict(final_lr_on_all_fof, newx = all_fof, type = "class")
first_step_fof_alpha = final_lr_on_all_fof$a0
first_step_fof_phage_beta = final_lr_on_all_fof$beta[1:phage_mut_cnt]
first_step_fof_bac_beta = final_lr_on_all_fof$beta[(phage_mut_cnt+1):(phage_mut_cnt+bac_mut_cnt)]

first_step_fof_phage_beta_l = data.frame(x = 1:phage_mut_cnt, y = 1, coef = first_step_fof_phage_beta)
first_step_fof_bac_beta_l = data.frame(x = 1:bac_mut_cnt, y = 1, coef = first_step_fof_bac_beta)

final_1st_fof_min = min(final_lr_on_all_fof$beta)
final_1st_fof_max = max(final_lr_on_all_fof$beta)
 
#step 1 logistic regression with sof predict(object,newdata,interval)
#final_1st_lr_on_all_sof_cv = cv.glmnet(all_sof, all_label, nfolds = 10, alpha = 1, family = "binomial", type.measure = "class", standardize = F)

# final_lr_on_all_sof = glmnet(all_sof, all_label, family = "binomial", alpha = 1, lambda = mean_1st_sof_lambda, standardize = F) # convergence error, lambda is too large?
final_lr_on_all_sof = glmnet(all_sof, all_label, family = "binomial", alpha = 1, lambda = 0, standardize = F) # convergence error, lambda is too large?

first_step_sof_yhat = predict(final_lr_on_all_sof, newx = all_sof, type = "class")
first_step_sof_alpha = final_lr_on_all_sof$a0
first_step_sof_beta_matrix = matrix(final_lr_on_all_sof$beta, nrow = bac_mut_cnt, ncol = phage_mut_cnt)
first_step_sof_beta_matrix_l = melt(first_step_sof_beta_matrix)
colnames(first_step_sof_beta_matrix_l) = c("bac_mut", "phage_mut", "coef")

final_1st_sof_min = min(final_lr_on_all_sof$beta)
final_1st_sof_max = max(final_lr_on_all_sof$beta)

#step 1 logistic regression with cof
#final_1st_lr_on_all_cof_cv = cv.glmnet(all_cof, all_label, nfolds = 10, alpha = 1, family = "binomial", type.measure = "class", standardize = F)
# final_lr_on_all_cof = glmnet(all_cof, all_label, family = "binomial", alpha = 1, lambda = mean_1st_cof_lambda, standardize = F)
final_lr_on_all_cof = glmnet(all_cof, all_label, family = "binomial", alpha = 1, lambda = 0, standardize = F) # convergence error, lambda is too large?

first_step_cof_yhat = predict(final_lr_on_all_cof, newx = all_cof, type = "class")
first_step_cof_alpha = final_lr_on_all_cof$a0
first_step_cof_phage_beta = final_lr_on_all_cof$beta[1:phage_mut_cnt]
first_step_cof_bac_beta = final_lr_on_all_cof$beta[(phage_mut_cnt+1):(phage_mut_cnt+bac_mut_cnt)]
first_step_cof_beta_matrix = matrix(final_lr_on_all_cof$beta[(phage_mut_cnt+bac_mut_cnt+1):length(final_lr_on_all_cof$beta)], nrow = bac_mut_cnt, ncol = phage_mut_cnt)

first_step_cof_phage_beta_l = data.frame(x = 1:phage_mut_cnt, y = 1, coef = first_step_cof_phage_beta)
first_step_cof_bac_beta_l = data.frame(x = 1:bac_mut_cnt, y = 1, coef = first_step_cof_bac_beta)
first_step_cof_beta_matrix_l = melt(first_step_cof_beta_matrix)
colnames(first_step_cof_beta_matrix_l) = c("bac_mut", "phage_mut", "coef")

final_1st_cof_min = min(final_lr_on_all_cof$beta)
final_1st_cof_max = max(final_lr_on_all_cof$beta)








#step2 linear wiht fof
# all pos data
pos_label = label_allgen[label_allgen > 0]
pos_label_log = log(pos_label)
pos_fof_bo = p_aa_b_aa_1st[label_allgen > 0, 178:196]
pos_fof_po = p_aa_b_aa_1st[label_allgen > 0, 1:177]
pos_fof = p_aa_b_aa_1st[label_allgen > 0, ]
pos_sof = p_aa_b_aa_2nd[label_allgen > 0, ]
pos_cof = cbind(pos_fof, pos_sof)

#step 2 linear regression with fof bac only
final_2nd_lm_on_pos_fof_bo = glmnet(pos_fof_bo, pos_label_log, family = "gaussian", alpha = 1, lambda = mean_2nd_fof_bo_lambda, standardize = F)

second_step_fof_bo_yhat  = predict(final_2nd_lm_on_pos_fof_bo, newx = pos_fof_bo)
second_step_fof_bo_yhat_rescale = exp(second_step_fof_bo_yhat)
second_step_fof_bo_alpha = final_2nd_lm_on_pos_fof_bo$a0
second_step_fof_bo_bac_beta = as.numeric(final_2nd_lm_on_pos_fof_bo$beta)

second_step_fof_bo_bac_beta_l = data.frame(x = 1:bac_mut_cnt, y = 1, coef = second_step_fof_bo_bac_beta)

final_2nd_fof_bo_min = min(final_2nd_lm_on_pos_fof_bo$beta)
final_2nd_fof_bo_max = max(final_2nd_lm_on_pos_fof_bo$beta)

#step 2 linear regression with fof phage only
final_2nd_lm_on_pos_fof_po = glmnet(pos_fof_po, pos_label_log, family = "gaussian", alpha = 1, lambda = mean_2nd_fof_po_lambda, standardize = F)

second_step_fof_po_yhat  = predict(final_2nd_lm_on_pos_fof_po, newx = pos_fof_po)
second_step_fof_po_yhat_rescale = exp(second_step_fof_po_yhat)
second_step_fof_po_alpha = final_2nd_lm_on_pos_fof_po$a0
second_step_fof_po_phage_beta = as.numeric(final_2nd_lm_on_pos_fof_po$beta)

second_step_fof_po_phage_beta_l = data.frame(x = 1:phage_mut_cnt, y = 1, coef = second_step_fof_po_phage_beta)

final_2nd_fof_po_min = min(final_2nd_lm_on_pos_fof_po$beta)
final_2nd_fof_po_max = max(final_2nd_lm_on_pos_fof_po$beta)


#step 2 linear regression with fof
final_2nd_lm_on_pos_fof = glmnet(pos_fof, pos_label_log, family = "gaussian", alpha = 1, lambda = mean_2nd_fof_lambda, standardize = F)

second_step_fof_yhat  = predict(final_2nd_lm_on_pos_fof, newx = pos_fof)
second_step_fof_yhat_rescale = exp(second_step_fof_yhat)
second_step_fof_alpha = final_2nd_lm_on_pos_fof$a0
second_step_fof_phage_beta = final_2nd_lm_on_pos_fof$beta[1:phage_mut_cnt]
second_step_fof_bac_beta = final_2nd_lm_on_pos_fof$beta[(phage_mut_cnt+1):(phage_mut_cnt+bac_mut_cnt)]

second_step_fof_phage_beta_l = data.frame(x = 1:phage_mut_cnt, y = 1, coef = second_step_fof_phage_beta)
second_step_fof_bac_beta_l = data.frame(x = 1:bac_mut_cnt, y = 1, coef = second_step_fof_bac_beta)

final_2nd_fof_min = min(final_2nd_lm_on_pos_fof$beta)
final_2nd_fof_max = max(final_2nd_lm_on_pos_fof$beta)

#step 2 linear regression with sof
final_2nd_lm_on_pos_sof = glmnet(pos_sof, pos_label, family = "gaussian", alpha = 1, lambda = mean_2nd_sof_lambda, standardize = F)

second_step_sof_yhat  = predict(final_2nd_lm_on_pos_sof, newx = pos_sof)
second_step_sof_yhat_rescale = exp(second_step_sof_yhat)
second_step_sof_alpha = final_2nd_lm_on_pos_sof$a0
second_step_sof_phage_beta = final_2nd_lm_on_pos_sof$beta[1:phage_mut_cnt]
second_step_sof_bac_beta = final_2nd_lm_on_pos_sof$beta[(phage_mut_cnt+1):(phage_mut_cnt+bac_mut_cnt)]
second_step_sof_beta_matrix = matrix(final_2nd_lm_on_pos_sof$beta, nrow = bac_mut_cnt, ncol = phage_mut_cnt)
second_step_sof_beta_matrix_l = melt(second_step_sof_beta_matrix)
colnames(second_step_sof_beta_matrix_l) = c("bac_mut", "phage_mut", "coef")

final_2nd_sof_min = min(final_2nd_lm_on_pos_sof$beta)
final_2nd_sof_max = max(final_2nd_lm_on_pos_sof$beta)

#step 2 linear regression with cof
final_2nd_lm_on_pos_cof = glmnet(pos_cof, pos_label, family = "gaussian", alpha = 1, lambda = mean_2nd_cof_lambda, standardize = F)

second_step_cof_yhat  = predict(final_2nd_lm_on_pos_cof, newx = pos_cof)
second_step_cof_yhat_rescale = exp(second_step_cof_yhat)
second_step_cof_alpha = final_2nd_lm_on_pos_cof$a0
second_step_cof_phage_beta = final_2nd_lm_on_pos_cof$beta[1:phage_mut_cnt]
second_step_cof_bac_beta = final_2nd_lm_on_pos_cof$beta[(phage_mut_cnt+1):(phage_mut_cnt+bac_mut_cnt)]
second_step_cof_beta_matrix = matrix(final_2nd_lm_on_pos_cof$beta[(phage_mut_cnt+bac_mut_cnt+1):length(final_lr_on_all_cof$beta)], nrow = bac_mut_cnt, ncol = phage_mut_cnt)

second_step_cof_phage_beta_l = data.frame(x = 1:phage_mut_cnt, y = 1, coef = second_step_cof_phage_beta)
second_step_cof_bac_beta_l = data.frame(x = 1:bac_mut_cnt, y = 1, coef = second_step_cof_bac_beta)
second_step_cof_beta_matrix_l = melt(second_step_cof_beta_matrix)
colnames(second_step_cof_beta_matrix_l) = c("bac_mut", "phage_mut", "coef")

final_2nd_cof_min = min(final_2nd_lm_on_pos_cof$beta)
final_2nd_cof_max = max(final_2nd_lm_on_pos_cof$beta)


#ses_GLMNET = ridge_se(pos_fof, pos_label, l_yhat, final_lm_on_pos_fof)

#barplot_df = data.frame(feature = factor(colnames(pos_fof), levels = colnames(pos_fof)), beta = as.numeric(final_lm_on_pos_fof$beta), se = ses_GLMNET)
# 
# pdf("p3_f4_final_lm_performance_pos.pdf", height = 4, width = 10, useDingbats = F)
# ggplot(barplot_df) + 
#   geom_point(aes(x = feature, y = beta), stat = "identity") +
#   geom_errorbar(aes(x=feature, ymin = beta - se, ymax = beta + se), width = 0.4) + 
#   geom_hline(yintercept = 0, linetype = 3, color = "red") + 
#   theme_bw() +
#   theme(panel.grid = element_blank(),
#         axis.text.x = element_text(angle=45,hjust=1))
# dev.off()

save.image("final_model_res3_FINAL.RData")

