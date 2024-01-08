mse <- function(truth, pred){
  mse = mean((truth - pred)^2)
  return(mse)
}

mae <- function(truth, pred){
  mae = mean(abs(truth-pred))
  return(mae)
}

mae_inflate <- function(truth, pred, inflate){
  mae_in = sum(abs(truth-pred))/(length(truth)+inflate)
  return(mae_in)
}

missC_error <- function(truth, pred){
  me = sum(truth != pred)/length(truth)
  return(me)
}

rf_model <- function(n, orig_label, pw_conc_f, pw_conv_f, pw_conc_p_conv_f, knowledge_sr_pair){
  
}

linear_model <- function(n, cut_level, orig_label, pw_conc_f, pw_conv_f, pw_conc_p_conv_f, knowledge_sr_pair){
  null_model_mae = NULL
  pw_conc_mae = NULL
  pw_conv_mae = NULL
  pw_conc_p_conv_mae = NULL
  
  cut_label = orig_label
  if(cut_level > 0){
    cut_label[cut_label > cut_level] = cut_level
  }
  
  training_idx = sample.int(length(cut_label), round(length(cut_label)*0.8), replace = T)
  test_idx = sample.int(length(cut_label), round(length(cut_label)*0.2), replace = T)
  
  training_label = cut_label[training_idx]
  test_label = cut_label[test_idx]
  
  pw_conc_training_feature = pw_conc_f[training_idx, ]
  pw_conc_test_feature = pw_conc_f[test_idx, ]
  
  pw_conv_training_feature = pw_conv_f[training_idx, ]
  pw_conv_test_feature = pw_conv_f[test_idx, ]
  
  pw_conc_p_conv_training_feature = pw_conc_p_conv_f[training_idx, ]
  pw_conc_p_conv_test_feature = pw_conc_p_conv_f[test_idx, ]
  
  training_mean = mean(cut_label[training_idx])
  
  training_idx_after_ignore = training_idx[which(!(training_idx %in% knowledge_sr_pair) == T)]
  training_knowledge_m_feature = pw_conc_p_conv_f[training_idx_after_ignore, ]
  training_knowledge_m_label = cut_label[training_idx_after_ignore]
  test_idx_ignore = which(test_idx %in% knowledge_sr_pair)
  test_knowledge_m_feature = pw_conc_p_conv_f[test_idx, ]
  
  # version 2 of knowledge, do not use merged feature
  training_knowledge_m_feature2 = pw_conv_f[training_idx_after_ignore, ]
  test_knowledge_m_feature2 = pw_conv_f[test_idx, ]
  
  ##########training phase
  cv_linear_pw_conc = cv.glmnet(pw_conc_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_pw_conv = cv.glmnet(pw_conv_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_pw_conc_p_conv = cv.glmnet(pw_conc_p_conv_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_knowledge_m_part = cv.glmnet(training_knowledge_m_feature, training_knowledge_m_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_knowledge_m_part2 = cv.glmnet(training_knowledge_m_feature2, training_knowledge_m_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  
  test_prediction_pw_conc = predict(cv_linear_pw_conc, pw_conc_test_feature, s = "lambda.min")
  test_prediction_pw_conv = predict(cv_linear_pw_conv, pw_conv_test_feature, s = "lambda.min")
  test_prediction_pw_conc_p_conv = predict(cv_linear_pw_conc_p_conv, pw_conc_p_conv_test_feature, s = "lambda.min")
  test_prediction_knowledge_m = predict(cv_linear_knowledge_m_part, test_knowledge_m_feature, s = "lambda.min")
  test_prediction_knowledge_m[test_idx_ignore] = 0
  test_prediction_knowledge_m2 = predict(cv_linear_knowledge_m_part2, test_knowledge_m_feature2, s = "lambda.min")
  test_prediction_knowledge_m2[test_idx_ignore] = 0
  
  test_ckneg_pw_conc = min(test_prediction_pw_conc)
  test_ckneg_pw_conv = min(test_prediction_pw_conv)
  test_ckneg_pw_conc_p_conv = min(test_prediction_pw_conc_p_conv)
  test_ckneg_pw_k1 = min(test_prediction_knowledge_m)
  test_ckneg_pw_k2 = min(test_prediction_knowledge_m2)
  
  test_mae_null = mae(test_label, rep(training_mean, length(test_label)))
  test_mae_pw_conc = mae(test_label, test_prediction_pw_conc)
  test_mae_pw_conv = mae(test_label, test_prediction_pw_conv)
  test_mae_pw_conc_p_conv = mae(test_label, test_prediction_pw_conc_p_conv)
  test_mae_knowledge_m = mae(test_label, test_prediction_knowledge_m)
  test_mae_knowledge_m2 = mae(test_label, test_prediction_knowledge_m2)
  
  return(c(null_mae = test_mae_null, pairwise_concatenation_mae = test_mae_pw_conc, pairwise_convolution_mae = test_mae_pw_conv, pairwise_conc_plus_conv_mae = test_mae_pw_conc_p_conv, knowledge_model_mae = test_mae_knowledge_m, knowledge_model_mae2 = test_mae_knowledge_m2, test_minimum_check_negative_conc = test_ckneg_pw_conc, test_minimum_check_negative_conv = test_ckneg_pw_conv, test_minimum_check_negative_conc_p_conv = test_ckneg_pw_conc_p_conv, test_minimum_check_negative_k1 = test_ckneg_pw_k1, test_minimum_check_negative_k2 = test_ckneg_pw_k2))
}


linear_model_test_hard_zero <- function(n, cut_level, orig_label, pw_conc_f, pw_conv_f, pw_conc_p_conv_f, knowledge_sr_pair){
  null_model_mae = NULL
  pw_conc_mae = NULL
  pw_conv_mae = NULL
  pw_conc_p_conv_mae = NULL
  
  cut_label = orig_label
  if(cut_level > 0){
    cut_label[cut_label > cut_level] = cut_level
  }
  
  training_idx = sample.int(length(cut_label), round(length(cut_label)*0.8), replace = T)
  test_idx = sample.int(length(cut_label), round(length(cut_label)*0.2), replace = T)
  
  training_label = cut_label[training_idx]
  test_label = cut_label[test_idx]
  
  pw_conc_training_feature = pw_conc_f[training_idx, ]
  pw_conc_test_feature = pw_conc_f[test_idx, ]
  
  pw_conv_training_feature = pw_conv_f[training_idx, ]
  pw_conv_test_feature = pw_conv_f[test_idx, ]
  
  pw_conc_p_conv_training_feature = pw_conc_p_conv_f[training_idx, ]
  pw_conc_p_conv_test_feature = pw_conc_p_conv_f[test_idx, ]
  
  training_mean = mean(cut_label[training_idx])
  
  training_idx_after_ignore = training_idx[which(!(training_idx %in% knowledge_sr_pair) == T)]
  training_knowledge_m_feature = pw_conc_p_conv_f[training_idx_after_ignore, ]
  training_knowledge_m_label = cut_label[training_idx_after_ignore]
  test_idx_ignore = which(test_idx %in% knowledge_sr_pair)
  test_knowledge_m_feature = pw_conc_p_conv_f[test_idx, ]
  
  # version 2 of knowledge, do not use merged feature
  training_knowledge_m_feature2 = pw_conv_f[training_idx_after_ignore, ]
  test_knowledge_m_feature2 = pw_conv_f[test_idx, ]
  
  ##########training phase
  cv_linear_pw_conc = cv.glmnet(pw_conc_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_pw_conv = cv.glmnet(pw_conv_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_pw_conc_p_conv = cv.glmnet(pw_conc_p_conv_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_knowledge_m_part = cv.glmnet(training_knowledge_m_feature, training_knowledge_m_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_knowledge_m_part2 = cv.glmnet(training_knowledge_m_feature2, training_knowledge_m_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  
  test_prediction_pw_conc = predict(cv_linear_pw_conc, pw_conc_test_feature, s = "lambda.min")
  test_prediction_pw_conc[test_prediction_pw_conc < 0] = 0
  test_prediction_pw_conv = predict(cv_linear_pw_conv, pw_conv_test_feature, s = "lambda.min")
  test_prediction_pw_conv[test_prediction_pw_conv < 0] = 0
  test_prediction_pw_conc_p_conv = predict(cv_linear_pw_conc_p_conv, pw_conc_p_conv_test_feature, s = "lambda.min")
  test_prediction_pw_conc_p_conv[test_prediction_pw_conc_p_conv < 0] = 0
  test_prediction_knowledge_m = predict(cv_linear_knowledge_m_part, test_knowledge_m_feature, s = "lambda.min")
  test_prediction_knowledge_m[test_prediction_knowledge_m < 0] = 0
  test_prediction_knowledge_m[test_idx_ignore] = 0
  test_prediction_knowledge_m2 = predict(cv_linear_knowledge_m_part2, test_knowledge_m_feature2, s = "lambda.min")
  test_prediction_knowledge_m2[test_prediction_knowledge_m2 < 0] = 0
  test_prediction_knowledge_m2[test_idx_ignore] = 0
  
  test_ckneg_pw_conc = min(test_prediction_pw_conc)
  test_ckneg_pw_conv = min(test_prediction_pw_conv)
  test_ckneg_pw_conc_p_conv = min(test_prediction_pw_conc_p_conv)
  test_ckneg_pw_k1 = min(test_prediction_knowledge_m)
  test_ckneg_pw_k2 = min(test_prediction_knowledge_m2)
  
  test_mae_null = mae(test_label, rep(training_mean, length(test_label)))
  test_mae_pw_conc = mae(test_label, test_prediction_pw_conc)
  test_mae_pw_conv = mae(test_label, test_prediction_pw_conv)
  test_mae_pw_conc_p_conv = mae(test_label, test_prediction_pw_conc_p_conv)
  test_mae_knowledge_m = mae(test_label, test_prediction_knowledge_m)
  test_mae_knowledge_m2 = mae(test_label, test_prediction_knowledge_m2)
  
  return(c(null_mae = test_mae_null, pairwise_concatenation_mae = test_mae_pw_conc, pairwise_convolution_mae = test_mae_pw_conv, pairwise_conc_plus_conv_mae = test_mae_pw_conc_p_conv, knowledge_model_mae = test_mae_knowledge_m, knowledge_model_mae2 = test_mae_knowledge_m2, test_minimum_check_negative_conc = test_ckneg_pw_conc, test_minimum_check_negative_conv = test_ckneg_pw_conv, test_minimum_check_negative_conc_p_conv = test_ckneg_pw_conc_p_conv, test_minimum_check_negative_k1 = test_ckneg_pw_k1, test_minimum_check_negative_k2 = test_ckneg_pw_k2))
}


linear_model_log_check <- function(n, cut_level, orig_label, pw_conc_f, pw_conv_f, pw_conc_p_conv_f, knowledge_sr_pair){
  null_model_mae = NULL
  pw_conc_mae = NULL
  pw_conv_mae = NULL
  pw_conc_p_conv_mae = NULL
  
  cut_label = log(orig_label)
  if(cut_level > 0){
    cut_label[cut_label > cut_level] = cut_level
  }
  
  training_idx = sample.int(length(cut_label), round(length(cut_label)*0.8), replace = T)
  test_idx = sample.int(length(cut_label), round(length(cut_label)*0.2), replace = T)
  
  training_label = cut_label[training_idx]
  test_label = cut_label[test_idx]
  
  pw_conc_training_feature = pw_conc_f[training_idx, ]
  pw_conc_test_feature = pw_conc_f[test_idx, ]
  
  pw_conv_training_feature = pw_conv_f[training_idx, ]
  pw_conv_test_feature = pw_conv_f[test_idx, ]
  
  pw_conc_p_conv_training_feature = pw_conc_p_conv_f[training_idx, ]
  pw_conc_p_conv_test_feature = pw_conc_p_conv_f[test_idx, ]
  
  training_mean = mean(cut_label[training_idx])
  
  training_idx_after_ignore = training_idx[which(!(training_idx %in% knowledge_sr_pair) == T)]
  training_knowledge_m_feature = pw_conc_p_conv_f[training_idx_after_ignore, ]
  training_knowledge_m_label = cut_label[training_idx_after_ignore]
  test_idx_ignore = which(test_idx %in% knowledge_sr_pair)
  test_knowledge_m_feature = pw_conc_p_conv_f[test_idx, ]
  
  # version 2 of knowledge, do not use merged feature
  training_knowledge_m_feature2 = pw_conv_f[training_idx_after_ignore, ]
  test_knowledge_m_feature2 = pw_conv_f[test_idx, ]
  
  ##########training phase
  cv_linear_pw_conc = cv.glmnet(pw_conc_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_pw_conv = cv.glmnet(pw_conv_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_pw_conc_p_conv = cv.glmnet(pw_conc_p_conv_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_knowledge_m_part = cv.glmnet(training_knowledge_m_feature, training_knowledge_m_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  cv_linear_knowledge_m_part2 = cv.glmnet(training_knowledge_m_feature2, training_knowledge_m_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae", standardize = F)
  
  
  test_prediction_pw_conc = exp(predict(cv_linear_pw_conc, pw_conc_test_feature, s = "lambda.min"))
  test_prediction_pw_conv = exp(predict(cv_linear_pw_conv, pw_conv_test_feature, s = "lambda.min"))
  test_prediction_pw_conc_p_conv = exp(predict(cv_linear_pw_conc_p_conv, pw_conc_p_conv_test_feature, s = "lambda.min"))
  test_prediction_knowledge_m = exp(predict(cv_linear_knowledge_m_part, test_knowledge_m_feature, s = "lambda.min"))
  test_prediction_knowledge_m[test_idx_ignore] = 0
  test_prediction_knowledge_m2 = exp(predict(cv_linear_knowledge_m_part2, test_knowledge_m_feature2, s = "lambda.min"))
  test_prediction_knowledge_m2[test_idx_ignore] = 0
  
  test_ckneg_pw_conc = min(test_prediction_pw_conc)
  test_ckneg_pw_conv = min(test_prediction_pw_conv)
  test_ckneg_pw_conc_p_conv = min(test_prediction_pw_conc_p_conv)
  test_ckneg_pw_k1 = min(test_prediction_knowledge_m)
  test_ckneg_pw_k2 = min(test_prediction_knowledge_m2)
  
  test_ckpos_pw_conc = max(test_prediction_pw_conc)
  test_ckpos_pw_conv = max(test_prediction_pw_conv)
  test_ckpos_pw_conc_p_conv = max(test_prediction_pw_conc_p_conv)
  test_ckpos_pw_k1 = max(test_prediction_knowledge_m)
  test_ckpos_pw_k2 = max(test_prediction_knowledge_m2)
  
  test_mae_null = mae(exp(test_label), rep(training_mean, length(test_label)))
  test_mae_pw_conc = mae(exp(test_label), test_prediction_pw_conc)
  test_mae_pw_conv = mae(exp(test_label), test_prediction_pw_conv)
  test_mae_pw_conc_p_conv = mae(exp(test_label), test_prediction_pw_conc_p_conv)
  test_mae_knowledge_m = mae(exp(test_label), test_prediction_knowledge_m)
  test_mae_knowledge_m2 = mae(exp(test_label), test_prediction_knowledge_m2)
  
  test_rmse_null = rmse(exp(test_label), rep(training_mean, length(test_label)))
  test_rmse_pw_conc = rmse(exp(test_label), test_prediction_pw_conc)
  test_rmse_pw_conv = rmse(exp(test_label), test_prediction_pw_conv)
  test_rmse_pw_conc_p_conv = rmse(exp(test_label), test_prediction_pw_conc_p_conv)
  test_rmse_knowledge_m = rmse(exp(test_label), test_prediction_knowledge_m)
  test_rmse_knowledge_m2 = rmse(exp(test_label), test_prediction_knowledge_m2)
  
  return(c(null_mae = test_mae_null, pairwise_concatenation_mae = test_mae_pw_conc, pairwise_convolution_mae = test_mae_pw_conv, pairwise_conc_plus_conv_mae = test_mae_pw_conc_p_conv, knowledge_model_mae = test_mae_knowledge_m, knowledge_model_mae2 = test_mae_knowledge_m2, pairwise_concatenation_rmse = test_rmse_pw_conc, pairwise_convolution_rmse = test_rmse_pw_conv, pairwise_conc_plus_conv_rmse = test_rmse_pw_conc_p_conv, knowledge_model_rmse = test_rmse_knowledge_m, knowledge_model_rmse2 = test_rmse_knowledge_m2, test_minimum_check_negative_conc = test_ckneg_pw_conc, test_minimum_check_negative_conv = test_ckneg_pw_conv, test_minimum_check_negative_conc_p_conv = test_ckneg_pw_conc_p_conv, test_minimum_check_negative_k1 = test_ckneg_pw_k1, test_minimum_check_negative_k2 = test_ckneg_pw_k2, test_maximum_check_positive_conc = test_ckpos_pw_conc, test_maximum_check_positive_conv = test_ckpos_pw_conv, test_maximum_check_positive_conc_p_conv = test_ckpos_pw_conc_p_conv, test_maximum_check_positive_k1 = test_ckpos_pw_k1, test_maximum_check_positive_k2 = test_ckpos_pw_k2))
}

linear_model_exp_link <- function(n, cut_level, orig_label, pw_conc_f, pw_conv_f, pw_conc_p_conv_f, knowledge_sr_pair){
  null_model_mae = NULL
  pw_conc_mae = NULL
  pw_conv_mae = NULL
  pw_conc_p_conv_mae = NULL
  
  cut_label = orig_label
  if(cut_level > 0){
    cut_label[cut_label > cut_level] = cut_level
  }
  
  training_idx = sample.int(length(cut_label), round(length(cut_label)*0.8), replace = T)
  test_idx = sample.int(length(cut_label), round(length(cut_label)*0.2), replace = T)
  
  training_label = cut_label[training_idx]
  test_label = cut_label[test_idx]
  
  pw_conc_training_feature = pw_conc_f[training_idx, ]
  pw_conc_test_feature = pw_conc_f[test_idx, ]
  
  pw_conv_training_feature = pw_conv_f[training_idx, ]
  pw_conv_test_feature = pw_conv_f[test_idx, ]
  
  pw_conc_p_conv_training_feature = pw_conc_p_conv_f[training_idx, ]
  pw_conc_p_conv_test_feature = pw_conc_p_conv_f[test_idx, ]
  
  training_mean = mean(cut_label[training_idx])
  
  training_idx_after_ignore = training_idx[which(!(training_idx %in% knowledge_sr_pair) == T)]
  training_knowledge_m_feature = pw_conc_p_conv_f[training_idx_after_ignore, ]
  training_knowledge_m_label = cut_label[training_idx_after_ignore]
  test_idx_ignore = which(test_idx %in% knowledge_sr_pair)
  test_knowledge_m_feature = pw_conc_p_conv_f[test_idx, ]
  
  # version 2 of knowledge, do not use merged feature
  training_knowledge_m_feature2 = pw_conv_f[training_idx_after_ignore, ]
  test_knowledge_m_feature2 = pw_conv_f[test_idx, ]
  
  ##########training phase
  cv_linear_pw_conc = cv.glmnet(pw_conc_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae")
  
  cv_linear_pw_conv = cv.glmnet(pw_conv_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae")
  
  cv_linear_pw_conc_p_conv = cv.glmnet(pw_conc_p_conv_training_feature, training_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae")
  
  cv_linear_knowledge_m_part = cv.glmnet(training_knowledge_m_feature, training_knowledge_m_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae")
  
  cv_linear_knowledge_m_part2 = cv.glmnet(training_knowledge_m_feature2, training_knowledge_m_label, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae")
  
  training_label_exp = log(training_label + 1e-8)
  cv_linear_pw_conc_p_conv_exp_link = cv.glmnet(pw_conc_p_conv_training_feature, training_label_exp, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae")
  training_knowledge_m_label_exp = log(training_knowledge_m_label + 1e-8)
  cv_linear_knowledge_m_part_exp_link = cv.glmnet(training_knowledge_m_feature, training_knowledge_m_label_exp, nfolds = 10, alpha = 0.5, family = "gaussian", type.measure = "mae")
  
  
  test_prediction_pw_conc = predict(cv_linear_pw_conc, pw_conc_test_feature, s = "lambda.min")
  test_prediction_pw_conv = predict(cv_linear_pw_conv, pw_conv_test_feature, s = "lambda.min")
  test_prediction_pw_conc_p_conv = predict(cv_linear_pw_conc_p_conv, pw_conc_p_conv_test_feature, s = "lambda.min")
  test_prediction_knowledge_m = predict(cv_linear_knowledge_m_part, test_knowledge_m_feature, s = "lambda.min")
  test_prediction_knowledge_m[test_idx_ignore] = 0
  test_prediction_knowledge_m2 = predict(cv_linear_knowledge_m_part2, test_knowledge_m_feature2, s = "lambda.min")
  test_prediction_knowledge_m2[test_idx_ignore] = 0
  test_prediction_pw_conc_p_conv_exp_link = predict(cv_linear_pw_conc_p_conv_exp_link, pw_conc_p_conv_test_feature, s = "lambda.min")
  test_prediction_pw_conc_p_conv_exp_link_recov = exp(test_prediction_pw_conc_p_conv_exp_link) - 1e-8
  test_prediction_knowledge_m_exp_link = predict(cv_linear_knowledge_m_part_exp_link, test_knowledge_m_feature, s = "lambda.min")
  test_prediction_knowledge_m_exp_link_recov = exp(test_prediction_knowledge_m_exp_link) - 1e-8
  test_prediction_knowledge_m_exp_link_recov[test_idx_ignore] = 0
  test_prediction_pw_conc_p_conv_reset = test_prediction_pw_conc_p_conv
  test_prediction_pw_conc_p_conv_reset[test_prediction_pw_conc_p_conv_reset < 0] = 0
  
  
  
  test_mae_null = mae(test_label, rep(training_mean, length(test_label)))
  test_mae_pw_conc = mae(test_label, test_prediction_pw_conc)
  test_mae_pw_conv = mae(test_label, test_prediction_pw_conv)
  test_mae_pw_conc_p_conv = mae(test_label, test_prediction_pw_conc_p_conv)
  test_mae_knowledge_m = mae(test_label, test_prediction_knowledge_m)
  test_mae_knowledge_m2 = mae(test_label, test_prediction_knowledge_m2)
  test_mae_pw_conc_p_conv_exp_link = mae(test_label, test_prediction_pw_conc_p_conv_exp_link_recov)
  test_mae_knowledge_m_exp_link = mae(test_label, test_prediction_knowledge_m_exp_link_recov)
  test_mae_pw_conc_p_conv_reset = mae(test_label, test_prediction_pw_conc_p_conv_reset)

  return(c(null_mae = test_mae_null, pairwise_concatenation_mae = test_mae_pw_conc, pairwise_convolution_mae = test_mae_pw_conv, pairwise_conc_plus_conv_mae = test_mae_pw_conc_p_conv, knowledge_model_mae = test_mae_knowledge_m, knowledge_model_mae2 = test_mae_knowledge_m2, test_mae_pw_conc_p_conv_exp_link = test_mae_pw_conc_p_conv_exp_link, test_mae_knowledge_m_exp_link = test_mae_knowledge_m_exp_link, test_mae_pw_conc_p_conv_reset = test_mae_pw_conc_p_conv_reset))
}


calc_aic <- function(fit){
  tLL = fit$nulldev - deviance(fit)
  k = fit$df
  n = fit$nobs
  AIC = -tLL + 2 * k
  return(AIC)
}

calc_aicc <- function(fit){
  tLL = fit$nulldev - deviance(fit)
  k = fit$df
  n = fit$nobs
  AICc = -tLL + 2 * k + 2 * k * (k + 1)/(n - k - 1)
  return(AICc)
}

calc_bic <- function(fit){
  tLL = fit$nulldev - deviance(fit)
  k = fit$df
  n = fit$nobs
  BIC = log(n) * k - tLL
  return(BIC)
}


bk_linear_model_simple <- function(n, cut_level, orig_label, feature, suppose_training){
  res_mae = NULL
  
  cut_label = orig_label
  if(cut_level > 0){
    cut_label[cut_label > cut_level] = cut_level
  }
  
  inflate_diff = suppose_training - round(length(cut_label)*0.8)
  
  training_idx = sample.int(length(cut_label), round(length(cut_label)*0.8), replace = T)
  test_idx = sample.int(length(cut_label), round(length(cut_label)*0.2), replace = T)
  
  training_label = cut_label[training_idx]
  test_label = cut_label[test_idx]
  
  training_feature = feature[training_idx, ]
  test_feature = feature[test_idx, ]
  
  
  ##########training phase
  cv_linear_model = cv.glmnet(training_feature, training_label, nfolds = 10, alpha = 1, family = "gaussian", type.measure = "mae")
  
  test_prediction = predict(cv_linear_model, test_feature, s = "lambda.min")
  
  test_mae = mae_inflate(test_label, test_prediction, inflate_diff)
  
  return(test_mae)
}


bk_knowledge_model <- function(bac_f, vir_f, labels, bac_exc_idx){
  bac_sample_exc_idx = which(bac_f[bac_exc_idx, ] == 1)
  bac_f_exc = bac_f[, -bac_sample_exc_idx]
  
  labels_exc = as.vector(labels[-bac_sample_exc_idx, ])
  
  suppose_training_size = round(length(labels)*0.8)
  
  conv_bac_vir = kronecker(t(vir_f), t(bac_f_exc))
  
  auto_right = length(bac_sample_exc_idx) * ncol(vir_f)
  
  cut1_res = sapply(1:200, linear_model_simple, 1, labels_exc, conv_bac_vir, suppose_training_size)
  cut2_res = sapply(1:200, linear_model_simple, 2, labels_exc, conv_bac_vir, suppose_training_size)
  cut5_res = sapply(1:200, linear_model_simple, 5, labels_exc, conv_bac_vir, suppose_training_size)
  cut_none_res = sapply(1:200, linear_model_simple, 0, labels_exc, conv_bac_vir, suppose_training_size)
  
  print(paste0("mean cut 1 is ", mean(cut1_res)))
  print(paste0("mean cut 2 is ", mean(cut2_res)))
  print(paste0("mean cut 5 is ", mean(cut5_res)))
  print(paste0("mean cut none is ", mean(cut_none_res)))
  
  return(rbind(cut1_res, cut2_res, cut5_res, cut_none_res))
}

get_initial_appear <- function(feature_matrix){
  init_appear_column = apply(feature_matrix, 1, function(x){min(which(x>0))})
  init_appear_names = colnames(feature_matrix)[init_appear_column]
  init_appear_days = as.numeric(gsub("_|-", "", substr(init_appear_names, 5, 6)))
  return(init_appear_day = init_appear_days)
}


get_linear_combine <- function(vec, m, ref_t){
  lmodel = lm(as.vector(vec) ~ 0 + m)
  pos_combine_list = which(lmodel$coefficients > 1e-4)
  neg_combine_list = which(lmodel$coefficients < -1e-4)
  pos_combine_list_anno = ref_t[pos_combine_list - 1, ]
  neg_combine_list_anno = ref_t[neg_combine_list - 1, ]
  combine_list = list(pos = pos_combine_list, neg = neg_combine_list, pos_anno = pos_combine_list_anno, neg_anno = neg_combine_list_anno, coeffs = lmodel$coefficients)
  #print(length(neg_combine_list))
  return(combine_list)
}


redist_coef_vector_vir <- function(weight_vec, mapping_list){
  redist_weights = rep(0, max(mapping_list$before))
  after_set = unique(mapping_list$after)
  for(i in after_set){
    sub_mapping = mapping_list[which(mapping_list$after == i), ]
    size_of_w = nrow(sub_mapping)
    avg_w = weight_vec[i] / size_of_w
    redist_weights[sub_mapping$before] = redist_weights[sub_mapping$before] + avg_w
  }
  return(redist_weights)
}


aggr_vir_wight_vec <- function(weight_vec){
  name_set = unique(names(weight_vec))
  new_vec = rep(0, length(name_set))
  names(new_vec) = name_set
  for(i in 1:length(name_set)){
    i_idx = which(names(weight_vec) == name_set[i])
    new_vec[i] = sum(weight_vec[i_idx])
  }
  return(new_vec)
}

redist_coef_weights <- function(weights, mapping_list){
  redist_weights = matrix(0, nrow = nrow(weights), ncol = max(mapping_list$before))
  after_set = unique(mapping_list$after)
  for(i in after_set){
    sub_mapping = mapping_list[which(mapping_list$after == i), ]
    size_of_w = nrow(sub_mapping)
    avg_w = weights[, i] / size_of_w
    redist_weights[, sub_mapping$before] = redist_weights[, sub_mapping$before] + avg_w
  }
  return(redist_weights)
}


get_dup_genome_by_col <- function(df){
  res_list = list()
  for(i in 1:ncol(df)){
    i_identical_col_idx_set = NULL
    for(j in 1:ncol(df)){
      if(identical(df[, i], df[, j])){
        i_identical_col_idx_set = c(i_identical_col_idx_set, j)
      }
    }
    if(length(i_identical_col_idx_set) == 1){
      i_identical_col_idx_set = NULL
    }
    res_list[[i]] = i_identical_col_idx_set
  }
  return(res_list)
}



pairwise_ttest_pval <- function(x, adj, addname){
  t_test_res = matrix(NA, nrow(x)-1, nrow(x)-1)
  for(i in 1:(nrow(x)-2)){
    for(j in (i+1):(nrow(x)-1)){
      cur_ttest = t.test(x[i, ], x[j, ])
      cur_ttest_p = cur_ttest$p.value
      if(adj){
        t_test_res[i, j] = ifelse(cur_ttest_p < 2.2e-16, "< 2.2e-16", cur_ttest_p)
        t_test_res[j, i] = ifelse(cur_ttest_p < 2.2e-16, "< 2.2e-16", cur_ttest_p)
      }
      else{
        t_test_res[i, j] = cur_ttest_p
        t_test_res[j, i] = cur_ttest_p
      }
    }
  }
  if(addname){
    colnames(t_test_res) = c("null", "model 1", "model 2", "model 3", "model 4")
    rownames(t_test_res) = c("null", "model 1", "model 2", "model 3", "model 4")
  }
  return(t_test_res)
}

# check_nonsyn <- function(anno_col){
#   tmp = strsplit(anno_col, "\\s")
#   tmp2 = sapply(tmp, function(x){
#     if(length(x) == 0){
#       return(1)
#     }
#     else if(length(x) >=2 & x[1] != "intergenic" & substr(x[1], 1, 1) != substr(x[1], nchar(x[1]), nchar(x[1]))){ # check for nonsynonymous mutations
#       return(1)
#     }
#     else{
#       return(0)
#     }
#   })
#   return(tmp2)
# }

check_nonsyn <- function(anno_col){
  tmp = strsplit(anno_col, "\\xa0")
  tmp2 = sapply(tmp, function(x){
    #z=unlist(strsplit(x, "\\xa0"))
    #  x=unlist(strsplit(x,'\\xa0'))[1]
    if(length(x) == 0){
      return(0)
    }
    else if(length(x) >=2 & x[1] != "intergenic" & x[1]!="coding" & substr(x[1], 1, 1) != substr(x[1], nchar(x[1]), nchar(x[1]))){ # check for nonsynonymous mutations
      return(1)
    }
    else{
      return(0)
    }
  })
  return(tmp2)
}

get_dup_genome_mapping <- function(uniq_cols, full_m){
  dup_res <- lapply(uniq_cols, function(uc){
    tmp = apply(full_m, 2, function(x){return(identical(x, full_m[, uc]))})
    if(sum(tmp) == 1){
      return()
    }
    else{
      return(names(tmp[which(tmp == T)]))
    }
  })
  
  dup_res <- dup_res[!sapply(dup_res, is.null)]
  return(dup_res)
}

check_genome_same_or_not <- function(g1, g2, sameset){
  any(sapply(sameset, function(x){
    g1 %in% x & g2 %in% x
  }))
}

get_phylo_cor_res <- function(phym, dup_mapping){
  phym[1, ] = phym[1, ] + .Machine$double.eps
  sample_names = colnames(phym)
  reslist = list()
  dup_cor = NULL
  nondup_cor = NULL
  cnt = 0
  for(i in 1:(length(sample_names)-1)){
    iphy = phym[, i]
    iname = sample_names[i]
    for(j in (i+1):length(sample_names)){
      jphy = phym[, j]
      jname = sample_names[j]
      ijcor = cor(iphy, jphy)
      if(check_genome_same_or_not(iname, jname, dup_mapping)){
        dup_cor = c(dup_cor, ijcor)
        if(ijcor < 0){
          print(paste0(iname, " and ", jname, " same genome but phenotype cor is ", ijcor))
        }
      }else{
        nondup_cor = c(nondup_cor, ijcor)
      }
      cnt = cnt + 1
    }
  }
  reslist[['dup_cor']] = dup_cor
  reslist[['nondup_cor']] = nondup_cor
  
  return(reslist)
}

get_phylo_cor_res_nosep <- function(phym){
  phym[1, ] = phym[1, ] + .Machine$double.eps
  sample_names = colnames(phym)
  res = NULL
  cnt = 0
  for(i in 1:(length(sample_names)-1)){
    iphy = phym[, i]
    iname = sample_names[i]
    for(j in (i+1):length(sample_names)){
      jphy = phym[, j]
      jname = sample_names[j]
      ijcor = cor(iphy, jphy)
      res = c(res, ijcor)
      cnt = cnt + 1
    }
  }
  
  return(res)
}

get_unique_phenom_by_avg <- function(phenom, bac_map, phage_map){
  ## avg row first
  for(i in 1:length(bac_map)){
    row_same_genome_name <- bac_map[[i]]
    row_same_genome_idx <- which(rownames(phenom) %in% row_same_genome_name)
    avg_pheno <- colMeans(phenom[row_same_genome_idx, ])
    phenom[row_same_genome_idx[1], ] <- avg_pheno
    phenom <- phenom[-row_same_genome_idx[2:length(row_same_genome_idx)], ]
  }
  
  ## avg column then
  for(j in 1:length(phage_map)){
    col_same_genome_name <- phage_map[[j]]
    col_same_genome_idx <- which(colnames(phenom) %in% col_same_genome_name)
    avg_pheno <- rowMeans(phenom[, col_same_genome_idx])
    phenom[, col_same_genome_idx[1]] <- avg_pheno
    phenom <- phenom[, -col_same_genome_idx[2:length(col_same_genome_idx)]]
  }
  
  return(phenom)
}

basemodel_pref <- function(all_labels){
  a = NULL
  for(i in seq(-0.1, 0.5, 0.01)){
    a1 = mae(all_labels, rep(i, length(all_labels)))
    a = c(a, a1)
  }
  plot(seq(-0.1, 0.5, 0.01), a, xlab = "predict all value to a specific number", ylab = "error")
  print(paste0("Minimum error is ", min(a)))
}

cv10_rf_model <- function(){
  foldidx = createFolds(all_labels, k = 10, list = T, returnTrain = T)
  cv10_res = sapply(1:10, cv_10_rf_core, foldidx)
}

cv_10_rf_core <- function(i, foldidx){
  train_idx = foldidx[[i]]
  test_idx = setdiff(1:length(all_labels), train_idx)
  
  train_idx_noknow = intersect(train_idx, knowledge_super_other)
  train_idx_know = intersect(train_idx, knowledge_super_resistant_pairs)
  test_idx_noknow = intersect(test_idx, knowledge_super_other)
  test_idx_know = intersect(test_idx, knowledge_super_resistant_pairs)
  
  ith_rf_model_1 = randomForest(x = pairwise_concatenation[train_idx, ], y = all_labels[train_idx])
  ith_rf_model_2 = randomForest(x = pairwise_convolution[train_idx, ], y = all_labels[train_idx])
  ith_rf_model_3 = randomForest(x = pairwise_conc_p_conv[train_idx, ], y = all_labels[train_idx])
  
  test_pred_res_m1 = predict(ith_rf_model_1, pairwise_concatenation[test_idx, ])
  test_pred_res_m2 = predict(ith_rf_model_2, pairwise_convolution[test_idx, ])
  test_pred_res_m3 = predict(ith_rf_model_3, pairwise_conc_p_conv[test_idx, ])
  
  mae1 = mae(all_labels[test_idx], test_pred_res_m1)
  mae2 = mae(all_labels[test_idx], test_pred_res_m2)
  mae3 = mae(all_labels[test_idx], test_pred_res_m3)
  
  rmse1 = rmse(all_labels[test_idx], test_pred_res_m1)
  rmse2 = rmse(all_labels[test_idx], test_pred_res_m2)
  rmse3 = rmse(all_labels[test_idx], test_pred_res_m3)
  
  training_error_mae = c(rf_1_train_mae = mae(ith_rf_model_1$predicted, ith_rf_model_1$y), rf_2_train_mae = mae(ith_rf_model_2$predicted, ith_rf_model_2$y), rf_3_train_mae = mae(ith_rf_model_3$predicted, ith_rf_model_3$y))
  training_error_rmse = c(rf_1_train_rmse = sqrt(ith_rf_model_1$mse[length(ith_rf_model_1$mse)]), rf_2_train_rmse = sqrt(ith_rf_model_2$mse[length(ith_rf_model_2$mse)]), rf_1p2_train_rmse = sqrt(ith_rf_model_3$mse[length(ith_rf_model_3$mse)]))
  
  test_error_mae = c(rf_1_test_mae = mae1, rf_2_test_mae = mae2, rf_1p2_test_mae = mae3)
  test_error_rmse = c(rf_1_test_rmse = rmse1, rf_2_test_rmse = rmse2, rf_3_test_rmse = rmse3)
  res = c(training_error_mae, training_error_rmse, test_error_mae, test_error_rmse)
  return(res)
}


bootstrap_rf_model <- function(x){
  train_idx = sample.int(length(all_labels), round(length(all_labels)*0.8), replace = T)
  test_idx = setdiff(1:length(all_labels), train_idx)
  
  train_idx_noknow = intersect(train_idx, knowledge_super_other)
  train_idx_know = intersect(train_idx, knowledge_super_resistant_pairs)
  test_idx_noknow = intersect(test_idx, knowledge_super_other)
  test_idx_know = intersect(test_idx, knowledge_super_resistant_pairs)
  
  
  ith_rf_model_1 = randomForest(x = pairwise_concatenation[train_idx, ], y = all_labels[train_idx])
  ith_rf_model_2 = randomForest(x = pairwise_convolution[train_idx, ], y = all_labels[train_idx])
  ith_rf_model_3 = randomForest(x = pairwise_conc_p_conv[train_idx, ], y = all_labels[train_idx])
  
  test_pred_res_m1 = predict(ith_rf_model_1, pairwise_concatenation[test_idx, ])
  test_pred_res_m2 = predict(ith_rf_model_2, pairwise_convolution[test_idx, ])
  test_pred_res_m3 = predict(ith_rf_model_3, pairwise_conc_p_conv[test_idx, ])
  
  mae1 = mae(all_labels[test_idx], test_pred_res_m1)
  mae2 = mae(all_labels[test_idx], test_pred_res_m2)
  mae3 = mae(all_labels[test_idx], test_pred_res_m3)
  
  rmse1 = rmse(all_labels[test_idx], test_pred_res_m1)
  rmse2 = rmse(all_labels[test_idx], test_pred_res_m2)
  rmse3 = rmse(all_labels[test_idx], test_pred_res_m3)
  
  training_error_mae = c(rf_1_train_mae = mae(ith_rf_model_1$predicted, ith_rf_model_1$y), rf_2_train_mae = mae(ith_rf_model_2$predicted, ith_rf_model_2$y), rf_3_train_mae = mae(ith_rf_model_3$predicted, ith_rf_model_3$y))
  training_error_rmse = c(rf_1_train_rmse = sqrt(ith_rf_model_1$mse[length(ith_rf_model_1$mse)]), rf_2_train_rmse = sqrt(ith_rf_model_2$mse[length(ith_rf_model_2$mse)]), rf_1p2_train_rmse = sqrt(ith_rf_model_3$mse[length(ith_rf_model_3$mse)]))
  
  test_error_mae = c(rf_1_test_mae = mae1, rf_2_test_mae = mae2, rf_1p2_test_mae = mae3)
  test_error_rmse = c(rf_1_test_rmse = rmse1, rf_2_test_rmse = rmse2, rf_3_test_rmse = rmse3)
  #res = rbind(training_error_rmse = training_error_rmse, test_error_mae = test_error_mae, test_error_rmse = test_error_rmse)
  res = c(training_error_mae, training_error_rmse, test_error_mae, test_error_rmse)
  
  return(res)
}

construct_eqtl_feature_table_1st <- function(pm, bm){
  firstorder = rbind(cbind(pm, matrix(0, nrow = nrow(pm), ncol = ncol(bm))), cbind(matrix(0, nrow = nrow(bm), ncol = ncol(pm)), bm))
  return(firstorder)
}

construct_eqtl_feature_table_1st_inworking <- function(pm, bm){
  #firstorder = rbind(cbind(pm, matrix(0, nrow = nrow(pm), ncol = ncol(bm))), cbind(matrix(0, nrow = nrow(bm), ncol = ncol(pm)), bm))
  firstorder = NULL
  rns = NULL
  for(i in 1:ncol(pm)){
    for(j in 1:ncol(bm)){
      rn = paste0(colnames(pm)[i], ":", colnames(bm)[j])
      rns = c(rns, rn)
      v = c(pm[, i], bm[, j])
      firstorder = rbind(firstorder, v)
    }
  }
  rownames(firstorder) = rns
  return(firstorder)
}

construct_eqtl_feature_table_2nd <- function(pm, bm){
  secondorder = kronecker(t(pm), t(bm), make.dimnames = T)
}

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

get_formular <- function(model){
  sprintf('y=%.4f+%s', coef(model)[1], paste(sprintf('%.4f*%s', coef(model)[-1], names(coef(model)[-1]) ), collapse ='+'))
}

perm_test <- function(tb){
  raw_rsqr <- summary(lm(tb[, 2] ~ tb[, 1]))$r.squared
  perm_rsqr <- NULL
  for(i in 1:500){
    iperm <- sample(tb[, 2], replace = F)
    iperm_rsqr <- summary(lm(iperm ~ tb[, 1]))$r.squared
    perm_rsqr <- c(perm_rsqr, iperm_rsqr)
  }
  q99 = quantile(perm_rsqr, 0.99)
  # hist(perm_rsqr, breaks = 20, col = "gray", main = "True R-square vs 500 permutated R-squares", xlab = "R-squared value", xlim = c(-1, 1), ylim = c(0, 100))
  # abline(v = raw_rsqr, col = "red")
  # abline(v = q99, col = "blue", lty = 2)
  # legend(-1.0, 100, c("True R-square", "99% quantile of permuted R-squares"), lty = c(1, 2), col = c("red", "blue"))
  # 
  
  
  perm_rsqr_df <- data.frame(a = perm_rsqr)
  
  plotobj <- ggplot(perm_rsqr_df, aes(x = a)) + 
    geom_histogram(bins = 25, fill = "gray") + 
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 100))  + 
    theme_bw() + 
    geom_vline(data = data.frame(xint = c(q99, raw_rsqr), color = factor(c( "99% quantile of permuted R-squares", "True R-square"))), aes(xintercept = xint, color = color), linetype="dashed", show.legend = T) +
    labs(x="R-square value", title="True R-square vs 500 permutated R-squares")
    # geom_vline(xintercept = q99, col = "blue", linetype="dashed", show.legend = T) + 
    # geom_vline(xintercept = raw_rsqr, col = "red", show.legend = T) + xlab("R-squared value") + 
    #ggtitle("True R-square vs 500 permutated R-squares") # + 
  # scale_color_hue(name = "", labels = c("99% quantile of permuted R-squares", "True R-square")) + 
  # theme(legend.box = 'horizontal', legend.position = c(-1, 60))
  return(plotobj)
}



runMatrixEqtl = function(snp_file_name, expression_file_name, pvOutputThreshold = 1e-2){
  require(MatrixEQTL)
  useModel = modelLINEAR
  SNP_file_name = snp_file_name
  expression_file_name = expression_file_name
  output_file_name = tempfile()
  pvOutputThreshold = pvOutputThreshold
  
  ## Load genotype data
  
  snps = SlicedData$new();
  snps$fileDelimiter = "\t";      # the TAB character
  snps$fileOmitCharacters = "NA"; # denote missing values;
  snps$fileSkipRows = 1;          # one row of column labels
  snps$fileSkipColumns = 1;       # one column of row labels
  snps$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  snps$LoadFile(SNP_file_name);
  
  ## Load gene expression data
  
  gene = SlicedData$new();
  gene$fileDelimiter = "\t";      # the TAB character
  gene$fileOmitCharacters = "NA"; # denote missing values;
  gene$fileSkipRows = 1;          # one row of column labels
  gene$fileSkipColumns = 1;       # one column of row labels
  gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
  gene$LoadFile(expression_file_name);
  
  me = Matrix_eQTL_engine(
    snps = snps,
    gene = gene,
    output_file_name = output_file_name,
    pvOutputThreshold = pvOutputThreshold,
    useModel = useModel, 
    verbose = TRUE,
    pvalue.hist = "qqplot",
    min.pv.by.genesnp = FALSE,
    noFDRsaveMemory = FALSE)
  
  return(me)
}

prepRunMatrixEqtl = function(gen_data, phen_data, outprefix){
  df1 = data.frame(id = rownames(gen_data), gen_data)
  if(is.null(rownames(phen_data))){
    df2 = data.frame(rbind(c("phenotype", phen_data)))
    colnames(df2) = colnames(df1)
  }else{
    df2 = data.frame(id = rownames(phen_data), phen_data)
  }
  print(dim(df1))
  print(dim(df2))
  gen_out = paste0("./eqtl_prep_folder/matrixeqtl_", outprefix, "_gen.txt")
  phen_out = paste0("./eqtl_prep_folder/matrixeqtl_", outprefix, "_phen.txt")
  write.table(df1, gen_out, sep = "\t", row.names = F, quote = F)
  write.table(df2, phen_out, sep = "\t", row.names = F, quote = F)
  p_me = runMatrixEqtl(snp_file_name = gen_out,
                       expression_file_name = phen_out,
                       pvOutputThreshold = 1)
}


summaryMatrixEqtl = function(me_res, title = ""){
  plot(me_res)
  mtext(title, side = 3, line = -2)
  return(me_res$all$eqtls)
}


addcol_2ndOrder_eqtlSummary_res <- function(summ_m){
  split_p_b = strsplit(summ_m[, 1], ":")
  p_mut = sapply(split_p_b, `[[`, 1)
  b_mut = sapply(split_p_b, `[[`, 2)
  added_m = cbind(summ_m, phage_mutation = p_mut, bac_mutation = b_mut)
  return(added_m)
}