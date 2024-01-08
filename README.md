# Genotype to phenotype inference model
This is a machine learning method to infer phage and bacterial mutations driving changes in infection phenotypes.

## Experimental data in /data
We use experimental data from Gupta et al., 2022 where E. coli B strain REL606 and phage λ strain cI26 were co-cultured for a 37-day period. Samples were taken on checkpoint days for sequencing and pairwise quantitative plaque assays as described in (Gupta A P. S., 2022). Table 1a and Table 1b contain genome wide changes of 50 bacterial host (descended from E. coli B strain REL606) and 44 phage (descended from λ strain cI26) strains including 18 and 176 unique mutations for the host and phage respectively. Table 2 contain the list of interactions of all phage-bacterial pairs including ancestors as the efficiency of a phage isolate infecting a derived host strain relative to that for infecting the ancestral strain (EOP values) yielding a 51 by 45 cross-infection matrix. 

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
