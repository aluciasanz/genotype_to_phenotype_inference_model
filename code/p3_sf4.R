options(stringsAsFactors = FALSE)
library(plyr)
library(ggplot2)
library(tidyr)
library(ggrepel)
PATH='./figures_11032022/'
# first step

#all_fof
#first_step_fof_phage_beta
#first_step_fof_bac_beta

bac_fof_mutname = rownames(bac_only_eqtl_m_allmut_allgen)
phage_fof_mutname = rownames(phage_only_eqtl_m_allmut_allgen)

### table mut name match

bac_mut_table = rbind(c('bac_an', NA, NA, NA, NA, "host ancestor indicator", NA), cbind(name = bac_fof_mutname[-1], bac_table[, c(1, 2, 54:56)]))
phage_mut_table = rbind(c('phage_an', NA, NA, NA, NA, "phage ancestor indicator", NA), cbind(name = phage_fof_mutname[-1], phage_table[, c(1, 2, 48:50)]))

all_mut_table = rbind(bac_mut_table, phage_mut_table)
all_mut_table$position = gsub(",", "", all_mut_table$position)

all_1st_fof_coef = rbind(data.frame(name = phage_fof_mutname[first_step_fof_phage_beta != 0], coef_val = first_step_fof_phage_beta[first_step_fof_phage_beta != 0], sample = "phage"), data.frame(name = bac_fof_mutname[first_step_fof_bac_beta != 0], coef_val = first_step_fof_bac_beta[first_step_fof_bac_beta != 0], sample = "host"))
all_1st_fof_coef_sort = all_1st_fof_coef[order(all_1st_fof_coef$coef_val), ]
all_1st_fof_coef_sort$x = 1:nrow(all_1st_fof_coef_sort)

# step 1 200 run min max table
step1_200_min_max_table = data.frame(name = rownames(all_logistic_res_allf_coef_fof_m), minV = apply(all_logistic_res_allf_coef_fof_m, 1, min), maxV = apply(all_logistic_res_allf_coef_fof_m, 1, max), low10 = apply(all_logistic_res_allf_coef_fof_m, 1, function(x) quantile(x, 0.1)), high10 = apply(all_logistic_res_allf_coef_fof_m, 1, function(x) quantile(x, 0.9)))

fof_1st_sort_imp_features = join(join(all_1st_fof_coef_sort, all_mut_table, by = "name", type = "left"), step1_200_min_max_table, by = "name", type = "left")

fof_1st_sort_imp_features_sub_same_sign_quantile = fof_1st_sort_imp_features[which(sign(fof_1st_sort_imp_features$coef_val) == sign(fof_1st_sort_imp_features$low10) & sign(fof_1st_sort_imp_features$coef_val) == sign(fof_1st_sort_imp_features$high10)), ]
#fof_1st_sort_imp_features_sub_same_sign_minmax = fof_1st_sort_imp_features[which(sign(fof_1st_sort_imp_features$coef_val) == sign(fof_1st_sort_imp_features$minV) & sign(fof_1st_sort_imp_features$coef_val) == sign(fof_1st_sort_imp_features$maxV)), ]


# fof_1st_sort_imp_features_processed = apply(fof_1st_sort_imp_features, 2, function(x){gsub(",", "", x)})
# fof_1st_sort_imp_features_processed = apply(fof_1st_sort_imp_features, 2, function(x){gsub("\\s*", "", x)})
# fof_1st_sort_imp_features_processed = apply(fof_1st_sort_imp_features_processed, 2, function(x){gsub("→|\\s*→\\s*", "->", x)})
# fof_1st_sort_imp_features_processed = apply(fof_1st_sort_imp_features_processed, 2, function(x){gsub("←|\\s*←\\s*", "<-", x)})

# write.table(fof_1st_sort_imp_features, file =paste0(PATH, "s2.step1_sort_nonzero_feature.tsv"), quote = F, sep = "\t", row.names = F)




step1_fof_coef_plot = ggplot(data = fof_1st_sort_imp_features, aes(x = x, y = coef_val, label= paste0(annotation,gene,sep="_"))) + 
  geom_point(aes(x = x, y = coef_val, shape = sample, color = sample), size = 2.0, stroke=1) + 
  scale_shape_manual(values = c(19, 5)) + 
  scale_color_manual(values = c("#91bfdb", "#fc8d59")) + 
  scale_x_continuous(breaks = c(1, 5, 10,15, 20, 25, 30,35,40,45,50,55,60,65,70)) + 
  scale_y_continuous(breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10))+
  geom_errorbar(aes(x = x, ymin = low10, ymax = high10, color = sample), width = 1) + 
  geom_hline(yintercept = 0, linetype = 3, color = "black") + 
  labs(x = "Rank", y = "Coefficient") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
    text = element_text(size = 20))
fof_1st_sort_imp_features_sig= rbind.data.frame(fof_1st_sort_imp_features[(fof_1st_sort_imp_features$high10< -0.0 & fof_1st_sort_imp_features$low10< -0.0),],fof_1st_sort_imp_features[(fof_1st_sort_imp_features$high10>0.0 & fof_1st_sort_imp_features$low10>0.0),])

s1=step1_fof_coef_plot+geom_text_repel(data = fof_1st_sort_imp_features[rownames(fof_1st_sort_imp_features_sig),],
                                    aes(label = paste0(annotation,gene,sep="")),
                                    size = 4,
                                    box.padding = unit(2, "lines"),
                                    point.padding = unit(0.01, "lines"),
                                    max.overlaps = 70) 


# ggsave(paste0(PATH,"supp_fig4_sig_1.svg"), plot = step1_fof_coef_plot,width =6,height = 5)
# ggsave(paste0("supp_fig4_sig_1.svg"), plot = s1, width =10,height = 5)

  # $Significant <- ifelse(fof_1st_sort_imp_features$minV > 0.001 | fof_1st_sort_imp_features$minV<-0.001, "Sig", "Not Sig")
set.seed(40)
step1_fof_coef_plot = ggplot(data = fof_1st_sort_imp_features, aes(x = x, y = coef_val,label= paste0(annotation,gene,sep="_"))) + 
  geom_point(aes(x = x, y = coef_val, shape = sample, color = sample), size = 2.0, stroke=1) + 
  scale_shape_manual(values = c(19, 5)) + 
  scale_color_manual(values = c("#91bfdb", "#fc8d59")) + 
  scale_x_continuous(breaks = c(1, 5, 10,15, 20, 25, 30,35,40,45,50,55,60)) + 
  scale_y_continuous(breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10))+
  geom_errorbar(aes(x = x, ymin = low10, ymax = high10, color = sample), width = 1) + 
  geom_hline(yintercept = 0, linetype = 3, color = "black") + 
  labs(x = "Rank", y = "Coefficient") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        text = element_text(size = 20))

annotations<-sub("\\xa0.*",'',fof_1st_sort_imp_features$annotation)
genes<-sub("\\xa0.*", "",fof_1st_sort_imp_features$gene)
fof_1st_sort_imp_features$final_annotation<- paste0(annotations,"\U2192" ,genes)

pt<-step1_fof_coef_plot+geom_text_repel(data = fof_1st_sort_imp_features[rownames(fof_1st_sort_imp_features_sig),],
    aes(label = final_annotation),
    size = 4,
    box.padding = unit(2, "lines"),
    point.padding = unit(0.01, "lines"),
    max.overlaps = 70
  )
ggsave(paste0("supp_fig4_sig_1.svg"), plot = pt,width =10,height = 5)
# second step

#second_step_fof_phage_beta
#second_step_fof_bac_beta

all_2nd_fof_coef = rbind(data.frame(name = phage_fof_mutname[second_step_fof_phage_beta != 0], coef_val = second_step_fof_phage_beta[second_step_fof_phage_beta != 0], sample = "phage"), data.frame(name = bac_fof_mutname[second_step_fof_bac_beta != 0], coef_val = second_step_fof_bac_beta[second_step_fof_bac_beta != 0], sample = "host"))
all_2nd_fof_coef_sort = all_2nd_fof_coef[order(all_2nd_fof_coef$coef_val), ]
all_2nd_fof_coef_sort$x = 1:nrow(all_2nd_fof_coef_sort)

step2_200_min_max_table = data.frame(name = rownames(boot_2nd_step_run_mae_res_all_coef_fof_m), minV = apply(boot_2nd_step_run_mae_res_all_coef_fof_m, 1, min), maxV = apply(boot_2nd_step_run_mae_res_all_coef_fof_m, 1, max), low10 = apply(boot_2nd_step_run_mae_res_all_coef_fof_m, 1, function(x) quantile(x, 0.1)), high10 = apply(boot_2nd_step_run_mae_res_all_coef_fof_m, 1, function(x) quantile(x, 0.9)))

fof_2nd_sort_imp_features = join(join(all_2nd_fof_coef_sort, all_mut_table, by = "name", type = "left"), step2_200_min_max_table, by = "name", type = "left")

fof_2nd_sort_imp_features_sub_same_sign_quantile = fof_2nd_sort_imp_features[which(sign(fof_2nd_sort_imp_features$coef_val) == sign(fof_2nd_sort_imp_features$low10) & sign(fof_2nd_sort_imp_features$coef_val) == sign(fof_2nd_sort_imp_features$high10)), ]
#fof_2nd_sort_imp_features_sub_same_sign_minmax = fof_2nd_sort_imp_features[which(sign(fof_2nd_sort_imp_features$coef_val) == sign(fof_2nd_sort_imp_features$minV) & sign(fof_2nd_sort_imp_features$coef_val) == sign(fof_2nd_sort_imp_features$maxV)), ]

# write.table(fof_2nd_sort_imp_features, file = paste0(PATH,"s3.step2_sort_nonzero_feature.tsv"), quote = F, sep = "\t", row.names = F)

fof_2nd_sort_imp_features_sig= rbind.data.frame(fof_2nd_sort_imp_features[(fof_2nd_sort_imp_features$high10< -0.0 & fof_2nd_sort_imp_features$low10< -0.0),],fof_2nd_sort_imp_features[(fof_2nd_sort_imp_features$high10>0 & fof_2nd_sort_imp_features$low10>0.0),])
# $Significant <- ifelse(fof_1st_sort_imp_features$minV > 0.001 | fof_1st_sort_imp_features$minV<-0.001, "Sig", "Not Sig")
set.seed(40)
step2_fof_coef_plot = ggplot(data = fof_2nd_sort_imp_features, aes(x = x, y = coef_val)) + 
  geom_point(aes(x = x, y = coef_val, shape = sample, color = sample), size = 2.0, stroke=1) + 
  scale_shape_manual(values = c(19, 5)) + 
  scale_color_manual(values = c("#91bfdb", "#fc8d59")) + 
  scale_x_continuous(breaks = c(1, 5, 10,15, 20, 25, 30,35,40,45,50,55,60,65,70)) + 
  scale_y_continuous(breaks=c(-10,-8,-6,-4,-2,0,2,4,6,8,10))+
  geom_errorbar(aes(x = x, ymin = low10, ymax = high10, color = sample), width = 1) + 
  geom_hline(yintercept = 0, linetype = 3, color = "black") + 
  labs(x = "Rank", y = "Coefficient") + 
  theme_bw() + 
  theme(panel.grid = element_blank(),
        text = element_text(size = 20))

annotations<-sub("\\xa0.*",'',fof_2nd_sort_imp_features$annotation)
genes<-sub("\\xa0.*", "",fof_2nd_sort_imp_features$gene)

fof_2nd_sort_imp_features$final_annotation<- paste0(annotations,"\U2192" ,genes)  

pt<-step2_fof_coef_plot +geom_text_repel(data = fof_2nd_sort_imp_features[rownames(fof_2nd_sort_imp_features_sig),],
                                          aes(label = final_annotation),
                                          size = 4,
                                          box.padding = unit(2, "lines"),
                                          point.padding = unit(0.01, "lines"),
                                          max.overlaps = 70
  )
# ggsave(paste0(PATH,"supp_fig4_sig_2bis.svg"), plot = pt,width =8,height = 5)
ggsave("supp_fig4_2.svg", plot = pt ,width =10, height = 5)

