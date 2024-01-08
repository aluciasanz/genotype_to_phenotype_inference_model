options(stringsAsFactors = F)

library(ggplot2)
library(stringr)
library(reshape2)

color_palette <- c("#bdbdbd", "gray0")
PATH='./figures_04252022/'

# raw_list =r read.csv("raw_infection_list_mod8-1.csv", header=T)

bac_order = rev(c(paste0("8-", 1:10), paste0("15-", 1:10), paste0("22-", 1:10), paste0("28-", 1:10), paste0("37-", 1:10)))

phg_order = c(paste0("8-", 1:11), paste0("15-", 1:11), paste0("22-", 1:11), paste0("28-", 1:11))

inter_row_idx = which((raw_list$Host != "A") & !str_detect(raw_list$Host, "^606") & (raw_list$Phage != "A"))

sub_list = raw_list[inter_row_idx, ]

fac_sub_list = sub_list

bac_name_factor = as.factor(bac_order)
phg_name_factor = as.factor(phg_order)

fac_sub_list$Host = factor(fac_sub_list$Host, bac_name_factor)
fac_sub_list$Phage = factor(fac_sub_list$Phage, phg_name_factor)

fac_sub_list$EOP_bin = as.factor(ifelse(fac_sub_list$EOP > 0, 1, 0))

sub_matrix = dcast(data = fac_sub_list, Host ~ Phage, value.var="EOP")
row.names(sub_matrix) = sub_matrix[, 1]
sub_matrix = sub_matrix[, -1]


full_plot = ggplot(fac_sub_list, aes(Phage, Host)) + geom_tile(aes(fill = EOP_bin), colour = "white") + scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0")) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.text.x = element_text(angle = 90))
full_plot

simp_plot = ggplot(fac_sub_list, aes(Phage, Host)) + geom_tile(aes(fill = EOP_bin), colour = "white") + scale_fill_manual(values = color_palette, name = "",  labels = c("0", ">0")) + theme_bw() +  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none', panel.border = element_rect(fill = NULL, color = 'white'))
simp_plot 

filename = "fig1_bw"

ggsave(paste0(PATH, filename, ".svg"),  width=7,height = 7, plot = simp_plot)
ggsave(paste0(PATH, filename, "_withname.svg"),  width=7,height = 7 , plot = full_plot)



