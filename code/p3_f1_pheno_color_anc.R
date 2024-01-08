options(stringsAsFactors = F)

library(ggplot2)
library(stringr)
library(reshape2)

color_palette <- c("#bdbdbd", "#fed976", "#feb24c", "#fd8d3c", "#fc4e2a", "#e31a1c", "#b10026")
PATH='./figures_04252022/'

# raw_list = read.csv("raw_infection_list_mod8-1.csv", header=T)

bac_order = rev(c("bac_an", paste0("8-", 1:10), paste0("15-", 1:10), paste0("22-", 1:10), paste0("28-", 1:10), paste0("37-", 1:10)))
bac_name = c("bac_an", paste0("8-", 1:10), paste0("15-", 1:10), paste0("22-", 1:10), paste0("28-", 1:10), paste0("37-", 1:10))
# 
phg_order = c("phage_an", paste0("8-", 1:11), paste0("15-", 1:11), paste0("22-", 1:11), paste0("28-", 1:11))
phage_name = c("phage_an", paste0("8-", 1:11), paste0("15-", 1:11), paste0("22-", 1:11), paste0("28-", 1:11))
# 
# inter_row_idx = which((raw_list$Host != "A") & !str_detect(raw_list$Host, "^606") & (raw_list$Phage != "A"))

sub_list = raw_list[inter_row_idx, c(1,2,6)]

list_add_ancestral = rbind(data.frame(Phage = phg_order, Host = "bac_an", EOP = 1), data.frame(Phage = "phage_an", Host = bac_order[1:(length(bac_order)-1)], EOP = 0))

fac_sub_list = rbind(sub_list, list_add_ancestral)

bac_name_factor = as.factor(bac_order)
phage_name_factor = as.factor(phg_order)

fac_sub_list$Host = factor(fac_sub_list$Host, bac_name_factor)
fac_sub_list$Phage = factor(fac_sub_list$Phage, phage_name_factor)

fac_sub_list$EOP_bin = as.factor(cut(fac_sub_list$EOP, c(0, .Machine$double.eps, 0.2, 0.4, 0.6, 0.8, 1, Inf), right = F))

full_plot = ggplot(fac_sub_list, aes(Phage, Host)) + geom_tile(aes(fill = EOP_bin), colour = "white") + scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1", "1 and +")) + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle = 45))
full_plot
truth_simp_plot = ggplot(fac_sub_list, aes(Phage, Host)) + geom_tile(aes(fill = EOP_bin), colour = "white") + scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0 - 0.2", "0.2 - 0.4", "0.4 - 0.6", "0.6 - 0.8", "0.8 - 1", "1 and +")) + theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none')
truth_simp_plot 

filename = "fig1_including_ancestor"

ggsave(paste0(PATH, filename, ".pdf"), plot = truth_simp_plot,width = 7, height = 7)
ggsave(paste0(PATH, filename, "_withname.pdf"), plot = full_plot, width = 7, height = 7)

# save the matrix

write.csv(fac_sub_list, paste0(PATH, filename,".csv"), row.names = F )
