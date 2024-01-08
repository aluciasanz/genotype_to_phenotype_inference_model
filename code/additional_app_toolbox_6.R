second_step_lm_pos_log <- function(n, label, fof, sof, cof, boot_n, err_measure){
  
  # bootstrap samples
  training_idx = sample.int(length(label), round(length(label)*boot_n), replace = T)
  test_idx = setdiff(1:length(label), unique(training_idx))
  
  trainY = label[training_idx]
  testY = label[test_idx]
  
  trainX_fof = fof[training_idx, ]
  testX_fof = fof[test_idx, ]
  
  trainX_sof = sof[training_idx, ]
  testX_sof = sof[test_idx, ]
  
  trainX_cof = cof[training_idx, ]
  testX_cof = cof[test_idx, ]
  
  trainY_mean = mean(trainY)
  
  # cross validation to select lambda
  cv_lm_fof = cv.glmnet(trainX_fof, trainY, nfolds = 10, alpha = 1, family = "gaussian", type.measure = err_measure, standardize = F)
  cv_lm_sof = cv.glmnet(trainX_sof, trainY, nfolds = 10, alpha = 1, family = "gaussian", type.measure = err_measure, standardize = F)
  cv_lm_cof = cv.glmnet(trainX_cof, trainY, nfolds = 10, alpha = 1, family = "gaussian", type.measure = err_measure, standardize = F)
  
  fof_lambda = cv_lm_fof$lambda.min
  sof_lambda = cv_lm_sof$lambda.min
  cof_lambda = cv_lm_cof$lambda.min
  
  # training
  # lm_fof = glmnet(trainX_fof, trainY, family = "gaussian", alpha = 0.5, lambda = 1, standardize = F)
  # lm_sof = glmnet(trainX_sof, trainY, family = "gaussian", alpha = 0.5, lambda = 1, standardize = F)
  # lm_cof = glmnet(trainX_cof, trainY, family = "gaussian", alpha = 0.5, lambda = 1, standardize = F)
  
  # prediction on training set
  train_pred_fof = predict(cv_lm_fof, trainX_fof, s = "lambda.min")
  train_pred_sof = predict(cv_lm_sof, trainX_sof, s = "lambda.min")
  train_pred_cof = predict(cv_lm_cof, trainX_cof, s = "lambda.min")
  
  # prediction on test set
  test_pred_fof = predict(cv_lm_fof, testX_fof, s = "lambda.min")
  test_pred_sof = predict(cv_lm_sof, testX_sof, s = "lambda.min")
  test_pred_cof = predict(cv_lm_cof, testX_cof, s = "lambda.min")
  
  # error calculation
  if(err_measure == "mae"){
    train_null_error = mae(trainY, rep(trainY_mean, length(trainY)))
    train_fof_error = mae(trainY, train_pred_fof)
    train_sof_error = mae(trainY, train_pred_sof)
    train_cof_error = mae(trainY, train_pred_cof)
    test_null_error = mae(testY, rep(trainY_mean, length(testY)))
    test_fof_error = mae(testY, test_pred_fof)
    test_sof_error = mae(testY, test_pred_sof)
    test_cof_error = mae(testY, test_pred_cof)
  }else{
    train_null_error = mse(trainY, rep(trainY_mean, length(trainY)))
    train_fof_error = mse(trainY, train_pred_fof)
    train_sof_error = mse(trainY, train_pred_sof)
    train_cof_error = mse(trainY, train_pred_cof)
    test_null_error = mse(testY, rep(trainY_mean, length(testY)))
    test_fof_error = mse(testY, test_pred_fof)
    test_sof_error = mse(testY, test_pred_sof)
    test_cof_error = mse(testY, test_pred_cof)
  }
  
  #return(c(train_null_error, train_fof_error, train_sof_error, train_cof_error, test_null_error, test_fof_error, test_sof_error, test_cof_error))
  return(c(train_null = train_null_error, train_fof = train_fof_error, train_sof = train_sof_error, train_cof = train_cof_error, test_null = test_null_error, test_fof = test_fof_error, test_sof = test_sof_error, test_cof = test_cof_error, fof_lambda = fof_lambda, sof_lambda = sof_lambda, cof_lambda = cof_lambda))
}

ridge_se <- function(xs,y,yhat,my_mod){
  # Note, you can't estimate an intercept here
  n <- dim(xs)[1]
  k <- dim(xs)[2]
  sigma_sq <- sum((y-yhat)^2)/ (n-k)
  lam <- my_mod$lambda
  #if(is.null(my_mod$lambda.min)==TRUE){lam <- 0}
  i_lams <- Matrix(diag(x=1,nrow=k,ncol=k),sparse=TRUE)
  xpx <- t(xs)%*%xs
  xpxinvplam <- solve(xpx+lam*i_lams)
  var_cov <- sigma_sq * (xpxinvplam %*% xpx %*% xpxinvplam)
  se_bs <- sqrt(diag(var_cov))
  print('NOTE: These standard errors are very biased.')
  return(se_bs)
}


first_step_logistic_all <- function(n, label, fof, sof, cof, boot_n, err_measure){
  training_idx = sample.int(length(label), round(length(label)*boot_n), replace = T)
  test_idx = setdiff(1:length(label), unique(training_idx))
  
  trainY = label[training_idx]
  testY = label[test_idx]
  
  trainX_fof = fof[training_idx, ]
  testX_fof = fof[test_idx, ]
  
  trainX_sof = sof[training_idx, ]
  testX_sof = sof[test_idx, ]
  
  trainX_cof = cbind(trainX_fof, trainX_sof)
  testX_cof = cbind(testX_fof, testX_sof)
  
  trainY_posrate_null = sum(label == 1)/length(label)
  
  # cross validation to select lambda
  cv_lm_fof = cv.glmnet(trainX_fof, trainY, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  cv_lm_sof = cv.glmnet(trainX_sof, trainY, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  cv_lm_cof = cv.glmnet(trainX_cof, trainY, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  
  fof_lambda = cv_lm_fof$lambda.min
  sof_lambda = cv_lm_sof$lambda.min
  cof_lambda = cv_lm_cof$lambda.min
  
  # training
  # lm_fof = glmnet(trainX_fof, trainY, family = "gaussian", alpha = 0.5, lambda = 1, standardize = F)
  # lm_sof = glmnet(trainX_sof, trainY, family = "gaussian", alpha = 0.5, lambda = 1, standardize = F)
  # lm_cof = glmnet(trainX_cof, trainY, family = "gaussian", alpha = 0.5, lambda = 1, standardize = F)
  
  # prediction on training set
  train_pred_fof = predict(cv_lm_fof, trainX_fof, s = "lambda.min", type = "class")
  train_pred_sof = predict(cv_lm_sof, trainX_sof, s = "lambda.min", type = "class")
  train_pred_cof = predict(cv_lm_cof, trainX_cof, s = "lambda.min", type = "class")
  
  # prediction on test set
  test_pred_fof = predict(cv_lm_fof, testX_fof, s = "lambda.min", type = "class")
  test_pred_sof = predict(cv_lm_sof, testX_sof, s = "lambda.min", type = "class")
  test_pred_cof = predict(cv_lm_cof, testX_cof, s = "lambda.min", type = "class")
  
  train_null_error = trainY_posrate_null
  train_fof_error = missC_error(trainY, train_pred_fof)
  train_sof_error = missC_error(trainY, train_pred_sof)
  train_cof_error = missC_error(trainY, train_pred_cof)
  test_null_error = missC_error(testY, ifelse(runif(length(testY)) > trainY_posrate_null, "0", "1"))
  test_fof_error = missC_error(testY, test_pred_fof)
  test_sof_error = missC_error(testY, test_pred_sof)
  test_cof_error = missC_error(testY, test_pred_cof)
  
  return(c(train_null = train_null_error, train_fof = train_fof_error, train_sof = train_sof_error, train_cof = train_cof_error, test_null = test_null_error, test_fof = test_fof_error, test_sof = test_sof_error, test_cof = test_cof_error, fof_lambda = fof_lambda, sof_lambda = sof_lambda, cof_lambda = cof_lambda))
}


first_step_logistic_all_sof_tune <- function(n, label, fof, sof, err_measure){
  boot_n = c(1, 1.2, 1.5)
  training_idx1 = sample.int(length(label), round(length(label)*boot_n[1]), replace = T)
  test_idx1 = setdiff(1:length(label), unique(training_idx1))
  
  training_idx2 = sample.int(length(label), round(length(label)*boot_n[2]), replace = T)
  test_idx2 = setdiff(1:length(label), unique(training_idx2))
  
  training_idx3 = sample.int(length(label), round(length(label)*boot_n[3]), replace = T)
  test_idx3 = setdiff(1:length(label), unique(training_idx3))
  
  trainY1 = label[training_idx1]
  testY1 = label[test_idx1]
  
  trainX_fof = fof[training_idx1, ]
  testX_fof = fof[test_idx1, ]
  
  trainX_sof1 = sof[training_idx1, ]
  testX_sof1 = sof[test_idx1, ]
  
  trainY2 = label[training_idx2]
  testY2 = label[test_idx2]
  
  trainX_sof2 = sof[training_idx2, ]
  testX_sof2 = sof[test_idx2, ]  
  
  trainY3 = label[training_idx3]
  testY3 = label[test_idx3]
  
  trainX_sof3 = sof[training_idx3, ]
  testX_sof3 = sof[test_idx3, ]
  
  # cross validation to select lambda
  cv_lm_fof = cv.glmnet(trainX_fof, trainY1, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  
  cv_lm_sof1 = cv.glmnet(trainX_sof1, trainY1, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  cv_lm_sof2 = cv.glmnet(trainX_sof2, trainY2, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  cv_lm_sof3 = cv.glmnet(trainX_sof3, trainY3, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  
  
  lm_sof1 = glmnet(trainX_sof1, trainY1, family = "binomial", alpha = 1, lambda = 0.005, standardize = F)
  lm_sof2 = glmnet(trainX_sof1, trainY1, family = "binomial", alpha = 1, lambda = 0.01, standardize = F)
  lm_sof3 = glmnet(trainX_sof1, trainY1, family = "binomial", alpha = 1, lambda = 0.02, standardize = F)
  

  # prediction on training set
  train_pred_fof = predict(cv_lm_fof, trainX_fof, s = "lambda.min", type = "class")
  train_pred_sof1 = predict(cv_lm_sof1, trainX_sof1, s = "lambda.min", type = "class")
  train_pred_sof2 = predict(cv_lm_sof2, trainX_sof2, s = "lambda.min", type = "class")
  train_pred_sof3 = predict(cv_lm_sof3, trainX_sof3, s = "lambda.min", type = "class")
  train_pred_sof4 = predict(lm_sof1, trainX_sof1, type = "class")
  train_pred_sof5 = predict(lm_sof2, trainX_sof1, type = "class")
  train_pred_sof6 = predict(lm_sof3, trainX_sof1, type = "class")
  

  # prediction on test set
  test_pred_fof = predict(cv_lm_fof, testX_fof, s = "lambda.min", type = "class")
  test_pred_sof1 = predict(cv_lm_sof1, testX_sof1, s = "lambda.min", type = "class")
  test_pred_sof2 = predict(cv_lm_sof2, testX_sof2, s = "lambda.min", type = "class")
  test_pred_sof3 = predict(cv_lm_sof3, testX_sof3, s = "lambda.min", type = "class")
  test_pred_sof4 = predict(lm_sof1, testX_sof1, type = "class")
  test_pred_sof5 = predict(lm_sof2, testX_sof1, type = "class")
  test_pred_sof6 = predict(lm_sof3, testX_sof1, type = "class")
  
  
  #error
  train_fof_error = missC_error(trainY1, train_pred_fof)
  train_sof_error1 = missC_error(trainY1, train_pred_sof1)
  train_sof_error2 = missC_error(trainY2, train_pred_sof2)
  train_sof_error3 = missC_error(trainY3, train_pred_sof3)
  train_sof_error4 = missC_error(trainY1, train_pred_sof4)
  train_sof_error5 = missC_error(trainY1, train_pred_sof5)
  train_sof_error6 = missC_error(trainY1, train_pred_sof6)
  
  test_fof_error = missC_error(testY1, test_pred_fof)
  test_sof_error1 = missC_error(testY1, test_pred_sof1)
  test_sof_error2 = missC_error(testY2, test_pred_sof2)
  test_sof_error3 = missC_error(testY3, test_pred_sof3)
  test_sof_error4 = missC_error(testY1, test_pred_sof4)
  test_sof_error5 = missC_error(testY1, test_pred_sof5)
  test_sof_error6 = missC_error(testY1, test_pred_sof6)
  

  return(c(train_fof_error, train_sof_error1, train_sof_error2, train_sof_error3, train_sof_error4, train_sof_error5, train_sof_error6, test_fof_error, test_sof_error1, test_sof_error2, test_sof_error3, test_sof_error4, test_sof_error5, test_sof_error6))
}



list_coef_bac_plot_n <- function(bac_coef_list, min_v, max_v){
  plotobj = ggplot(bac_coef_list, aes(x=x, y=y)) + 
    geom_tile(aes(x=x, y=y, fill = coef), col = "white") + 
    #xlab("Bac Mutations") + 
    ylab("Weights") +
    scale_fill_gradientn(colours = c("darkred", "red", "white","green", "darkgreen"), values = scales::rescale(c(min_v, max_v)), guide = "colorbar", limits=c(min_v, max_v), breaks=linspace(min_v, max_v, 5)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank())
}


first_step_logistic_all_5model <- function(n, label, fof_b, fof_p, fof, sof, cof, boot_n, err_measure){
  training_idx = sample.int(length(label), round(length(label)*boot_n), replace = T)
  test_idx = setdiff(1:length(label), unique(training_idx))
  
  test_idx_bined = rep(0, length(label))
  test_idx_bined[test_idx] = 1
  
  trainY = label[training_idx]
  testY = label[test_idx]
  
  trainX_fof_bo = fof_b[training_idx, ]
  testX_fof_bo = fof_b[test_idx, ]
  
  trainX_fof_po = fof_p[training_idx, ]
  testX_fof_po = fof_p[test_idx, ]
  
  trainX_fof = fof[training_idx, ]
  testX_fof = fof[test_idx, ]
  
  trainX_sof = sof[training_idx, ]
  testX_sof = sof[test_idx, ]
  
  trainX_cof = cbind(trainX_fof, trainX_sof)
  testX_cof = cbind(testX_fof, testX_sof)
  
  trainY_posrate_null = sum(label == 1)/length(label)
  
  # cross validation to select lambda
  cv_lm_fof_bo = cv.glmnet(trainX_fof_bo, trainY, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  cv_lm_fof_po = cv.glmnet(trainX_fof_po, trainY, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  cv_lm_fof = cv.glmnet(trainX_fof, trainY, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  cv_lm_sof = cv.glmnet(trainX_sof, trainY, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  cv_lm_cof = cv.glmnet(trainX_cof, trainY, nfolds = 10, alpha = 1, family = "binomial", type.measure = err_measure, standardize = F)
  
  fof_bo_lambda = cv_lm_fof_bo$lambda.min
  fof_po_lambda = cv_lm_fof_po$lambda.min
  fof_lambda = cv_lm_fof$lambda.min
  sof_lambda = cv_lm_sof$lambda.min
  cof_lambda = cv_lm_cof$lambda.min
  
  # get coeff
  lm_fof_bo = glmnet(trainX_fof_bo, trainY, family = "binomial", alpha = 1, lambda = cv_lm_fof_bo$lambda.min, standardize = F)
  lm_fof_po = glmnet(trainX_fof_po, trainY, family = "binomial", alpha = 1, lambda = cv_lm_fof_po$lambda.min, standardize = F)
  lm_fof = glmnet(trainX_fof, trainY, family = "binomial", alpha = 1, lambda = cv_lm_fof$lambda.min, standardize = F)
  lm_sof = glmnet(trainX_sof, trainY, family = "binomial", alpha = 1, lambda = cv_lm_sof$lambda.min, standardize = F)
  lm_cof = glmnet(trainX_cof, trainY, family = "binomial", alpha = 1, lambda = cv_lm_cof$lambda.min, standardize = F)

  coef_fof_bo = lm_fof_bo$beta
  coef_fof_po = lm_fof_po$beta
  coef_fof = lm_fof$beta
  coef_sof = lm_sof$beta
  coef_cof = lm_cof$beta
  
  # prediction on training set
  train_pred_fof_bo = predict(cv_lm_fof_bo, trainX_fof_bo, s = "lambda.min", type = "class")
  train_pred_fof_po = predict(cv_lm_fof_po, trainX_fof_po, s = "lambda.min", type = "class")
  train_pred_fof = predict(cv_lm_fof, trainX_fof, s = "lambda.min", type = "class")
  train_pred_sof = predict(cv_lm_sof, trainX_sof, s = "lambda.min", type = "class")
  train_pred_cof = predict(cv_lm_cof, trainX_cof, s = "lambda.min", type = "class")
  
  # prediction on test set
  test_pred_fof_bo = predict(cv_lm_fof_bo, testX_fof_bo, s = "lambda.min", type = "class")
  test_pred_fof_po = predict(cv_lm_fof_po, testX_fof_po, s = "lambda.min", type = "class")
  test_pred_fof = predict(cv_lm_fof, testX_fof, s = "lambda.min", type = "class")
  test_pred_sof = predict(cv_lm_sof, testX_sof, s = "lambda.min", type = "class")
  test_pred_cof = predict(cv_lm_cof, testX_cof, s = "lambda.min", type = "class")
  
  train_null_error = trainY_posrate_null
  train_fof_bo_error = missC_error(trainY, train_pred_fof_bo)
  train_fof_po_error = missC_error(trainY, train_pred_fof_po)
  train_fof_error = missC_error(trainY, train_pred_fof)
  train_sof_error = missC_error(trainY, train_pred_sof)
  train_cof_error = missC_error(trainY, train_pred_cof)
  test_null_error = missC_error(testY, ifelse(runif(length(testY)) > trainY_posrate_null, "0", "1"))
  test_fof_bo_error = missC_error(testY, test_pred_fof_bo)
  test_fof_po_error = missC_error(testY, test_pred_fof_po)
  test_fof_error = missC_error(testY, test_pred_fof)
  test_sof_error = missC_error(testY, test_pred_sof)
  test_cof_error = missC_error(testY, test_pred_cof)
  
  test_pred_result = cbind(test_pred_fof_bo, test_pred_fof_po, test_pred_fof, test_pred_sof, test_pred_cof)
  colnames(test_pred_result) = c('fof_bo','fof_po','fof','sof','cof')
  
  coef_trained = list(coef_fof_bo = coef_fof_bo, coef_fof_po = coef_fof_po, coef_fof = coef_fof, coef_sof = coef_sof, coef_cof = coef_cof)
  
  stats = c(train_null = train_null_error, train_fof_bo = train_fof_bo_error, train_fof_po = train_fof_po_error, train_fof = train_fof_error, train_sof = train_sof_error, train_cof = train_cof_error, test_null = test_null_error, test_fof_bo = test_fof_bo_error, test_fof_po = test_fof_po_error, test_fof = test_fof_error, test_sof = test_sof_error, test_cof = test_cof_error, fof_bo_lambda = fof_bo_lambda, fof_po_lambda = fof_po_lambda, fof_lambda = fof_lambda, sof_lambda = sof_lambda, cof_lambda = cof_lambda)
  
  results = list(stats = stats, test_idx = test_idx, test_idx_bined = test_idx_bined, test_pred_result = test_pred_result, coef_trained = coef_trained)
  
  return(results)
}


second_step_lm_pos_log_5model <- function(n, label, fof_b, fof_p, fof, sof, cof, boot_n, err_measure){
  # bootstrap samples
  training_idx = sample.int(length(label), round(length(label)*boot_n), replace = T)
  test_idx = setdiff(1:length(label), unique(training_idx))
  
  test_idx_bined = rep(0, length(label))
  test_idx_bined[test_idx] = 1
  
  trainY = label[training_idx]
  testY = label[test_idx]
  
  trainX_fof_bo = fof_b[training_idx, ]
  testX_fof_bo = fof_b[test_idx, ]
  
  trainX_fof_po = fof_p[training_idx, ]
  testX_fof_po = fof_p[test_idx, ]
  
  trainX_fof = fof[training_idx, ]
  testX_fof = fof[test_idx, ]
  
  trainX_sof = sof[training_idx, ]
  testX_sof = sof[test_idx, ]
  
  trainX_cof = cof[training_idx, ]
  testX_cof = cof[test_idx, ]
  
  trainY_mean = mean(trainY)
  
  # cross validation to select lambda
  cv_lm_fof_bo = cv.glmnet(trainX_fof_bo, trainY, nfolds = 10, alpha = 1, family = "gaussian", type.measure = err_measure, standardize = F)
  cv_lm_fof_po = cv.glmnet(trainX_fof_po, trainY, nfolds = 10, alpha = 1, family = "gaussian", type.measure = err_measure, standardize = F)
  cv_lm_fof = cv.glmnet(trainX_fof, trainY, nfolds = 10, alpha = 1, family = "gaussian", type.measure = err_measure, standardize = F)
  cv_lm_sof = cv.glmnet(trainX_sof, trainY, nfolds = 10, alpha = 1, family = "gaussian", type.measure = err_measure, standardize = F)
  cv_lm_cof = cv.glmnet(trainX_cof, trainY, nfolds = 10, alpha = 1, family = "gaussian", type.measure = err_measure, standardize = F)
  
  fof_bo_lambda = cv_lm_fof_bo$lambda.min
  fof_po_lambda = cv_lm_fof_po$lambda.min
  fof_lambda = cv_lm_fof$lambda.min
  sof_lambda = cv_lm_sof$lambda.min
  cof_lambda = cv_lm_cof$lambda.min
  
  # get training coef
  lm_fof_bo = glmnet(trainX_fof, trainY, family = "gaussian", alpha = 1, lambda = fof_bo_lambda, standardize = F)
  lm_fof_po = glmnet(trainX_fof, trainY, family = "gaussian", alpha = 1, lambda = fof_po_lambda, standardize = F)
  lm_fof = glmnet(trainX_fof, trainY, family = "gaussian", alpha = 1, lambda = fof_lambda, standardize = F)
  lm_sof = glmnet(trainX_sof, trainY, family = "gaussian", alpha = 1, lambda = sof_lambda, standardize = F)
  lm_cof = glmnet(trainX_cof, trainY, family = "gaussian", alpha = 1, lambda = cof_lambda, standardize = F)
  
  coef_fof_bo = lm_fof_bo$beta
  coef_fof_po = lm_fof_po$beta
  coef_fof = lm_fof$beta
  coef_sof = lm_sof$beta
  coef_cof = lm_cof$beta
  
  # prediction on training set
  train_pred_fof_bo = predict(cv_lm_fof_bo, trainX_fof_bo, s = "lambda.min")
  train_pred_fof_po = predict(cv_lm_fof_po, trainX_fof_po, s = "lambda.min")
  train_pred_fof = predict(cv_lm_fof, trainX_fof, s = "lambda.min")
  train_pred_sof = predict(cv_lm_sof, trainX_sof, s = "lambda.min")
  train_pred_cof = predict(cv_lm_cof, trainX_cof, s = "lambda.min")
  
  # prediction on test set
  test_pred_fof_bo = predict(cv_lm_fof_bo, testX_fof_bo, s = "lambda.min")
  test_pred_fof_po = predict(cv_lm_fof_po, testX_fof_po, s = "lambda.min")
  test_pred_fof = predict(cv_lm_fof, testX_fof, s = "lambda.min")
  test_pred_sof = predict(cv_lm_sof, testX_sof, s = "lambda.min")
  test_pred_cof = predict(cv_lm_cof, testX_cof, s = "lambda.min")
  
  # error calculation
  if(err_measure == "mae"){
    train_null_error = mae(trainY, rep(trainY_mean, length(trainY)))
    train_fof_bo_error = mae(trainY, train_pred_fof_bo)
    train_fof_po_error = mae(trainY, train_pred_fof_po)
    train_fof_error = mae(trainY, train_pred_fof)
    train_sof_error = mae(trainY, train_pred_sof)
    train_cof_error = mae(trainY, train_pred_cof)
    test_null_error = mae(testY, rep(trainY_mean, length(testY)))
    test_fof_bo_error = mae(testY, test_pred_fof_bo)
    test_fof_po_error = mae(testY, test_pred_fof_po)
    test_fof_error = mae(testY, test_pred_fof)
    test_sof_error = mae(testY, test_pred_sof)
    test_cof_error = mae(testY, test_pred_cof)
  }else{
    train_null_error = mse(trainY, rep(trainY_mean, length(trainY)))
    train_fof_bo_error = mse(trainY, train_pred_fof_bo)
    train_fof_po_error = mse(trainY, train_pred_fof_po)
    train_fof_error = mse(trainY, train_pred_fof)
    train_sof_error = mse(trainY, train_pred_sof)
    train_cof_error = mse(trainY, train_pred_cof)
    test_null_error = mse(testY, rep(trainY_mean, length(testY)))
    test_fof_bo_error = mse(testY, test_pred_fof_bo)
    test_fof_po_error = mse(testY, test_pred_fof_po)
    test_fof_error = mse(testY, test_pred_fof)
    test_sof_error = mse(testY, test_pred_sof)
    test_cof_error = mse(testY, test_pred_cof)
  }
  
  test_pred_result = cbind(test_pred_fof_bo, test_pred_fof_po, test_pred_fof, test_pred_sof, test_pred_cof)
  colnames(test_pred_result) = c('fof_bo','fof_po','fof','sof','cof')
  
  coef_trained = list(coef_fof_bo = coef_fof_bo, coef_fof_po = coef_fof_po, coef_fof = coef_fof, coef_sof = coef_sof, coef_cof = coef_cof)
  
  stats = c(train_null = train_null_error, train_fof_bo = train_fof_bo_error, train_fof_po = train_fof_po_error, train_fof = train_fof_error, train_sof = train_sof_error, train_cof = train_cof_error, test_null = test_null_error, test_fof_bo = test_fof_bo_error, test_fof_po = test_fof_po_error, test_fof = test_fof_error, test_sof = test_sof_error, test_cof = test_cof_error, fof_bo_lambda = fof_bo_lambda, fof_po_lambda = fof_po_lambda, fof_lambda = fof_lambda, sof_lambda = sof_lambda, cof_lambda = cof_lambda)
  
  results = list(stats = stats, test_idx = test_idx, test_idx_bined = test_idx_bined, test_pred_result = test_pred_result, coef_trained = coef_trained)
  
  return(results)
}