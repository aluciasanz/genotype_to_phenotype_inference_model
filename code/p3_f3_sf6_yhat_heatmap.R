options(stringsAsFactors = F)
#setwd("~/Documents/2019coevolution_geno_pheno/code_and_data_paper3/figs_als")

library(ggplot2)
library(stringr)
library(reshape2)

color_palette <- c("#bdbdbd", "gray0")
PATH='./figures_04252022/'
#comment this if you have it already in the environment
#raw_list = read.csv("raw_infection_list_mod8-1.csv", header=T)

bac_order = rev(c("bac_an", paste0("8-", 1:10), paste0("15-", 1:10), paste0("22-", 1:10), paste0("28-", 1:10), paste0("37-", 1:10)))

phg_order = c("phage_an", paste0("8-", 1:11), paste0("15-", 1:11), paste0("22-", 1:11), paste0("28-", 1:11))

inter_row_idx = which((raw_list$Host != "A") & !str_detect(raw_list$Host, "^606") & (raw_list$Phage != "A"))

sub_list = raw_list[inter_row_idx, c(1,2,6)]

list_add_ancestral = rbind(data.frame(Phage = phg_order, Host = "bac_an", EOP = 1), data.frame(Phage = "phage_an", Host = bac_order[1:(length(bac_order)-1)], EOP = 0))

fac_sub_list = rbind(sub_list, list_add_ancestral)

bac_name_factor = as.factor(bac_order)
phg_name_factor = as.factor(phg_order)

fac_sub_list$Host = factor(fac_sub_list$Host, bac_name_factor)
fac_sub_list$Phage = factor(fac_sub_list$Phage, phg_name_factor)

fac_sub_list$EOP_bin = as.factor(ifelse(fac_sub_list$EOP > 0, 1, 0))

raw_bin_pheno_m = ifelse(p_allgen_b_allgen_pheno_m > 0, 1, 0)

full_plot = ggplot(fac_sub_list, aes(Phage, Host)) + geom_tile(aes(fill = EOP_bin), colour = "white") + scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle = 90))
full_plot
truth_simp_plot = ggplot(fac_sub_list, aes(Phage, Host)) + geom_tile(aes(fill = EOP_bin), colour = "white") + scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0")) + theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none')
truth_simp_plot

filename = "fig3_truth_y"

ggsave(paste0(filename, ".svg"),width = 7, height = 7, plot = truth_simp_plot)
ggsave(paste0(filename, "_withname.svg"),width = 7, height = 7, plot = full_plot)


bac_name = c("bac_an", paste0("8-", 1:10), paste0("15-", 1:10), paste0("22-", 1:10), paste0("28-", 1:10), paste0("37-", 1:10))
phage_name = c("phage_an", paste0("8-", 1:11), paste0("15-", 1:11), paste0("22-", 1:11), paste0("28-", 1:11))

bac_order = rev(bac_name)
phage_order = phage_name

bac_name_factor = as.factor(bac_order)
phage_name_factor = as.factor(phage_order)

first_step_yhat_heatmap <- function(pred_yhat){
  first_step_yhat_m = matrix(as.numeric(as.character(pred_yhat)), nrow = bac_cnt, ncol = phage_cnt)
  rownames(first_step_yhat_m) = bac_name
  colnames(first_step_yhat_m) = phage_name
  first_step_yhat_m_l = melt(first_step_yhat_m)
  colnames(first_step_yhat_m_l) = c("bac", "phage", "val")
  first_step_yhat_m_l$val_bin = as.factor(ifelse(first_step_yhat_m_l$val == 0, 0, 1))
  
  first_step_yhat_m_l$bac = factor(first_step_yhat_m_l$bac, bac_name_factor)
  first_step_yhat_m_l$phage = factor(first_step_yhat_m_l$phage, phage_name_factor)
  
  # full_plot = ggplot(first_step_yhat_m_l, aes(phage, bac)) + 
  #   geom_tile(aes(fill = val_bin), colour = "white") + 
  #   scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0")) + 
  #   theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle = 90))
  
  simp_plot = ggplot(first_step_yhat_m_l, aes(phage, bac)) + 
    geom_tile(aes(fill = val_bin), colour = "white") + 
    scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0")) + 
    theme_bw() +  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none')
  
  return(simp_plot)
}

# first step yhats
first_step_fof_bo_yhat_plot = first_step_yhat_heatmap(first_step_fof_bo_yhat)
first_step_fof_bo_yhat_plot 

first_step_fof_po_yhat_plot = first_step_yhat_heatmap(first_step_fof_po_yhat)
first_step_fof_po_yhat_plot

first_step_fof_yhat_plot = first_step_yhat_heatmap(first_step_fof_yhat)
first_step_fof_yhat_plot

first_step_sof_yhat_plot = first_step_yhat_heatmap(first_step_sof_yhat)
first_step_sof_yhat_plot

first_step_cof_yhat_plot = first_step_yhat_heatmap(first_step_cof_yhat)
first_step_cof_yhat_plot

ggsave(paste0(PATH,"supp_fig6_fof_bo_yhat.svg"), plot = first_step_fof_bo_yhat_plot,width = 7, height = 7)
ggsave(paste0(PATH,"supp_fig6_fof_po_yhat.svg"), plot = first_step_fof_po_yhat_plot,width = 7, height = 7)
ggsave(paste0(PATH,"fig3_fof_yhat.svg"), plot = first_step_fof_yhat_plot,width = 7, height = 7)
ggsave(paste0(PATH,"fig3_sof_yhat.svg"), plot = first_step_sof_yhat_plot,width = 7, height = 7)
ggsave(paste0(PATH,"fig3_cof_yhat.svg"), plot = first_step_cof_yhat_plot,width = 7, height = 7)

