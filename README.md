# Genotype to phenotype inference model
## How to run the machine learning framework
### step 1. 
```
run approximation_uniq_nonsyn.R
```
### step 2. 
```
run all_feature_setup.R
```
### step 3.
```
run eqtl_trial.R
```
### step 4.(WARNING: code may take ~4h)
```
run first_step_logistic_regression.R
```
or 
```
load 'first_step_all_feature_logistic_run_res3.RData' in R. 
```
### step 5. (WARNING: code may take ~4h)
```
run second_step_linear_model.R
```
or
```
load 'second_step_all_feature_run_res3.RData'  
```
### step 6.
```
run final_model_using_all_data.R
```
or
```
load 'final_model_res3.RData' in R.
```
## Plot the results of the model.
### Figure 1: 
```
run 	p3_f1_mut_heatmap.R, 	p3_f1_pheno_blackwhite.R, 	p3_f1_pheno_color_anc.R, 	p3_f1_pheno_color.R 
```
### Figure 2:
```
run p3_f2_sf2.R 
```
### Figure 3, 4: 
```
run	p3_f3_sf6_yhat_heatmap.R, 	p3_f4_sf57_yhat_heatmap.R, p3_f34_sf567_coef_heatmap_plot.R
```
### Figure 5:
```
run p3_sf4.R
```
### S1,3,Figure:
```
run p3_sf13.R and stitch
```
