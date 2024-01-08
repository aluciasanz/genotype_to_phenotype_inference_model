options(stringsAsFactors = F)
#setwd("~/Documents/2019coevolution_geno_pheno/code_and_data_paper3/figs_als")
PATH='./figures_04252022/'
library(ggplot2)
library(stringr)
library(reshape2)

color_palette <- c("white", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")
color_palette_raw <- c("white", "#ccebc5", "#7bccc4", "#4eb3d3", "#2b8cbe", "#0868ac", "#084081")
#comment this if you have it already in the environment
#raw_list = read.csv("raw_infection_list_mod8-1.csv", header=T)

bac_order = rev(c("bac_an", paste0("8-", 1:10), paste0("15-", 1:10), paste0("22-", 1:10), paste0("28-", 1:10), paste0("37-", 1:10)))
bac_name = c("bac_an", paste0("8-", 1:10), paste0("15-", 1:10), paste0("22-", 1:10), paste0("28-", 1:10), paste0("37-", 1:10))

phg_order = c("phage_an", paste0("8-", 1:11), paste0("15-", 1:11), paste0("22-", 1:11), paste0("28-", 1:11))
phage_name = c("phage_an", paste0("8-", 1:11), paste0("15-", 1:11), paste0("22-", 1:11), paste0("28-", 1:11))

inter_row_idx = which((raw_list$Host != "A") & !str_detect(raw_list$Host, "^606") & (raw_list$Phage != "A"))

sub_list = raw_list[inter_row_idx, c(1,2,6)]

list_add_ancestral = rbind(data.frame(Phage = phg_order, Host = "bac_an", EOP = 1), data.frame(Phage = "phage_an", Host = bac_order[1:(length(bac_order)-1)], EOP = 0))

fac_sub_list = rbind(sub_list, list_add_ancestral)

bac_name_factor = as.factor(bac_order)
phage_name_factor = as.factor(phg_order)

fac_sub_list$Host = factor(fac_sub_list$Host, bac_name_factor)
fac_sub_list$Phage = factor(fac_sub_list$Phage, phage_name_factor)

fac_sub_list_log = fac_sub_list

fac_sub_list$EOP_bin = as.factor(cut(fac_sub_list$EOP, c(0, .Machine$double.eps, 0.2, 0.4, 0.6, 0.8, 1, Inf), right = F))

fac_sub_list_log$EOP[fac_sub_list_log$EOP == 0] = -9999
fac_sub_list_log$EOP[fac_sub_list_log$EOP > 0] = log(fac_sub_list_log$EOP[fac_sub_list_log$EOP > 0])
fac_sub_list_log$EOP_bin = as.factor(cut(fac_sub_list_log$EOP, c(-Inf, -9, -6.5, -4, -1.5, 1, 3.5, Inf), right = F))

full_plot = ggplot(fac_sub_list, aes(Phage, Host)) + geom_tile(aes(fill = EOP_bin), colour = "white") + scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1", "1 and +")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle = 45))
full_plot
truth_simp_plot = ggplot(fac_sub_list, aes(Phage, Host)) + geom_tile(aes(fill = EOP_bin), colour = "white") + scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1", "1 and +")) + theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none')
truth_simp_plot
full_plot_log = ggplot(fac_sub_list_log, aes(Phage, Host)) + geom_tile(aes(fill = EOP_bin), colour = "white") + scale_fill_manual(values = color_palette_raw, name = "",  labels = c("NA", "[-9, -6.5)", "[-6.5, -4)", "[-4, -1.5)", "[-1.5, 1)", "[1, 3.5)", "3.5+")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle = 45))
full_plot_log
truth_simp_plot_log = ggplot(fac_sub_list_log, aes(Phage, Host)) + geom_tile(aes(fill = EOP_bin), colour = "white") + scale_fill_manual(values = color_palette_raw, name = "",  labels = c("NA", "[-9, -6.5)", "[-6.5, -4)", "[-4, -1.5)", "[-1.5, 1)", "[1, 3.5)", "3.5+")) + theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none')
truth_simp_plot_log
filename = "fig4_truth_y"
ggsave(paste0(PATH,filename, ".svg"), plot = truth_simp_plot,width = 7, height = 7)
ggsave(paste0(PATH,filename, "_withname.svg"), plot = full_plot,width = 7, height = 7)

filename = "supp_fig5_truth_y_log"
ggsave(paste0(PATH,filename, ".svg"), plot = truth_simp_plot_log,width = 7, height = 7)
ggsave(paste0(PATH,filename, "_withname.svg"), plot = full_plot_log,width = 7, height = 7)


###################################################
#### plots for fof_bo, fof_po, fof, sof, cof
###################################################
template_m = acast(fac_sub_list, Host~Phage, value.var="EOP")
template_m = apply(template_m,2,rev)

second_step_yhat_heatmap <- function(second_step_yhat_rescale){
  second_step_yhat_rescale_m = template_m
  rownames(second_step_yhat_rescale_m) = bac_name
  colnames(second_step_yhat_rescale_m) = phage_name
  second_step_yhat_rescale_m[second_step_yhat_rescale_m > 0] = second_step_yhat_rescale
  
  second_step_yhat_rescale_m_l = melt(second_step_yhat_rescale_m)
  colnames(second_step_yhat_rescale_m_l) = c("bac", "phage", "val")
  second_step_yhat_rescale_m_l$val_bin = as.factor(cut(second_step_yhat_rescale_m_l$val, c(0, .Machine$double.eps, 0.2, 0.4, 0.6, 0.8, 1, Inf), right = F))
  second_step_yhat_rescale_m_l$col = cut(second_step_yhat_rescale_m_l$val, c(0, .Machine$double.eps, 0.2, 0.4, 0.6, 0.8, 1, Inf), right = F, labels = color_palette)
  
  
  second_step_yhat_rescale_m_l$bac = factor(second_step_yhat_rescale_m_l$bac, bac_name_factor)
  second_step_yhat_rescale_m_l$phage = factor(second_step_yhat_rescale_m_l$phage, phage_name_factor)
  
  # full_plot = ggplot(second_step_yhat_rescale_m_l, aes(phage, bac)) + geom_tile(aes(fill = val_bin), colour = "white") +  scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1", "1 and +")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle = 45))
  
  simp_plot = ggplot(second_step_yhat_rescale_m_l, aes(phage, bac)) + 
    geom_tile(aes(fill = val_bin), colour = "white") + 
    scale_fill_manual(values = levels(factor(second_step_yhat_rescale_m_l$col))) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none')
  return(simp_plot)
}

second_step_yhat_log_heatmap <- function(second_step_yhat){
  second_step_yhat_log_m = template_m
  rownames(second_step_yhat_log_m) = bac_name
  colnames(second_step_yhat_log_m) = phage_name
  second_step_yhat_log_m[second_step_yhat_log_m == 0] = -9999
  second_step_yhat_log_m[second_step_yhat_log_m > 0] = second_step_yhat
  
  second_step_yhat_log_m_l = melt(second_step_yhat_log_m)
  colnames(second_step_yhat_log_m_l) = c("bac", "phage", "val")
  second_step_yhat_log_m_l$val_bin = as.factor(cut(second_step_yhat_log_m_l$val, c(-Inf, -9, -6.5, -4, -1.5, 1, 3.5, Inf), right = F))
  second_step_yhat_log_m_l$col = cut(second_step_yhat_log_m_l$val, c(-Inf, -9, -6.5, -4, -1.5, 1, 3.5, Inf), right = F, labels = color_palette_raw)
  
  second_step_yhat_log_m_l$bac = factor(second_step_yhat_log_m_l$bac, bac_name_factor)
  second_step_yhat_log_m_l$phage = factor(second_step_yhat_log_m_l$phage, phage_name_factor)
  
  simp_plot = ggplot(second_step_yhat_log_m_l, aes(phage, bac)) + 
    geom_tile(aes(fill = val_bin), colour = "white") + 
    scale_fill_manual(values = levels(factor(second_step_yhat_log_m_l$col))) +
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none')
  return(simp_plot)
}


# second step yhats
second_step_fof_bo_yhat_log_plot = second_step_yhat_log_heatmap(second_step_fof_bo_yhat)
second_step_fof_bo_yhat_plot = second_step_yhat_heatmap(second_step_fof_bo_yhat_rescale)

second_step_fof_po_yhat_log_plot = second_step_yhat_log_heatmap(second_step_fof_po_yhat)
second_step_fof_po_yhat_plot = second_step_yhat_heatmap(second_step_fof_po_yhat_rescale)

second_step_fof_yhat_log_plot = second_step_yhat_log_heatmap(second_step_fof_yhat)
second_step_fof_yhat_plot = second_step_yhat_heatmap(second_step_fof_yhat_rescale)

second_step_sof_yhat_log_plot = second_step_yhat_log_heatmap(second_step_sof_yhat)
second_step_sof_yhat_plot = second_step_yhat_heatmap(second_step_sof_yhat_rescale)

second_step_cof_yhat_log_plot = second_step_yhat_log_heatmap(second_step_cof_yhat)
second_step_cof_yhat_plot = second_step_yhat_heatmap(second_step_cof_yhat_rescale)


ggsave(paste0(PATH,"supp_fig7_fof_bo_yhat_log.pdf"), plot = second_step_fof_bo_yhat_log_plot, width = 7, height = 7)
ggsave(paste0(PATH,"supp_fig7_fof_bo_yhat.pdf"), plot = second_step_fof_bo_yhat_plot, width = 7, height = 7)

ggsave(paste0(PATH,"supp_fig7_fof_po_yhat_log.pdf"), plot = second_step_fof_po_yhat_log_plot, width = 7, height = 7)
ggsave(paste0(PATH,"supp_fig7_fof_po_yhat.pdf"), plot = second_step_fof_po_yhat_plot, width = 7, height = 7)

ggsave(paste0(PATH,"supp_fig5_fof_yhat_log.pdf"), plot = second_step_fof_yhat_log_plot, width = 7, height = 7)
ggsave(paste0(PATH,"fig4_fof_yhat.pdf"), plot = second_step_fof_yhat_plot, width = 7, height = 7)

ggsave(paste0(PATH,"supp_fig5_sof_yhat_log.pdf"), plot = second_step_sof_yhat_log_plot, width = 7, height = 7)
ggsave(paste0(PATH,"fig4_sof_yhat.pdf"), plot = second_step_sof_yhat_plot, width = 7, height = 7)

ggsave(paste0(PATH,"supp_fig5_cof_yhat_log.pdf"), plot = second_step_cof_yhat_log_plot, width = 7, height = 7)
ggsave(paste0(PATH,"fig4_cof_yhat.pdf"), plot = second_step_cof_yhat_plot, width = 7, height = 7)
