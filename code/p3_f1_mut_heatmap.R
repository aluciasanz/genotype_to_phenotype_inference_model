options(stringsAsFactors = F)
library(stringr)
library(reshape2)
library(XML)
library(ggplot2)
library(svglite) # sudo apt-get install r-cran-svglite
library(tidyverse)
library(devtools)
#########################################
########### Figure 1 P1 #################
#########################################

color_palette <- c("gray90", "gray10")
PATH='./figures_04252022/'

################ bac mut plot
# bac_table_doc = htmlParse("bac_anno_50.html", encoding = "UTF-8")
# bac_table = data.frame(readHTMLTable(bac_table_doc, stringsAsFactors = FALSE, check.names = F))
# colnames(bac_table) = gsub("-", "_", bac_table[1, ])
# bac_table = bac_table[which(bac_table[, 1] != "position"), ]
# rownames(bac_table) = 1:nrow(bac_table)
# bac_feature_raw = bac_table[, 3:52]
# bac_feature_binary = ifelse(bac_feature_raw == "" | bac_feature_raw == "?", 0, 1)

bac_feature_binary_l = melt(bac_feature_binary)
bac_feature_binary_l$value = factor(bac_feature_binary_l$value, levels = c(0, 1))

bac_order = rev(c(paste0("B_D_8_", 1:10), paste0("B_D_15_", 1:10), paste0("B_D_22_", 1:10), paste0("B_D_28_", 1:10), paste0("B_D_37_", 1:10)))

bac_name_factor = as.factor(bac_order)
bac_feature_binary_l$Var2 = factor(bac_feature_binary_l$Var2, bac_name_factor)

bac_plot = ggplot(bac_feature_binary_l, aes(Var1, Var2)) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_manual(values = color_palette, name = "",  labels = c("contain", "not contain")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none', axis.line.x = element_blank(), axis.line.y = element_blank())
bac_plot 
ggsave(paste0(PATH,"fig1_hm_bac_mut.svg"), width=7,height = 7 , plot = bac_plot)


bac_feature_binary_add_an_l = melt(bac_only_eqtl_m_allmut_allgen)
bac_feature_binary_add_an_l$value = factor(bac_feature_binary_add_an_l$value, levels = c(0, 1))

bac_order_an = rev(c("an_b", paste0("B_D_8_", 1:10), paste0("B_D_15_", 1:10), paste0("B_D_22_", 1:10), paste0("B_D_28_", 1:10), paste0("B_D_37_", 1:10)))
bac_name_an_factor = as.factor(bac_order_an)
bac_feature_binary_add_an_l$Var2 = factor(bac_feature_binary_add_an_l$Var2, bac_name_an_factor)

bac_plot_an = ggplot(bac_feature_binary_add_an_l, aes(Var1, Var2)) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_manual(values = color_palette, name = "",  labels = c("contain", "not contain")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none', axis.line.x = element_blank(), axis.line.y = element_blank())

bac_plot_an 

ggsave(paste0(PATH,"fig1_hm_bac_mut_anc.svg"),  width=7,height = 7, plot = bac_plot_an)

############### phage mut plot
# phage_table_doc = htmlParse("phage_all_annotated_cut50.html", encoding = "UTF-8")
# phage_table = data.frame(readHTMLTable(phage_table_doc, stringsAsFactors = FALSE, check_names = F))
# colnames(phage_table) = gsub("-", "_", phage_table[1, ])
# phage_table = phage_table[which(phage_table[, 1] != "position"), ]
# rownames(phage_table) = 1:nrow(phage_table)
# phage_feature_raw = phage_table[, 3:46]
# phage_feature_binary = ifelse(phage_feature_raw == "" | phage_feature_raw == "?", 0, 1)

phage_feature_binary_l = melt(phage_feature_binary)
phage_feature_binary_l$value = factor(phage_feature_binary_l$value, levels = c(0, 1))

row_order = rev(sort(unique(phage_feature_binary_l$Var1)))
row_order_factor = as.factor(row_order)
phage_feature_binary_l$Var1 = factor(phage_feature_binary_l$Var1, row_order_factor)

phg_order = c(paste0("P_D_8_", 1:11), paste0("P_D_15_", 1:11), paste0("P_D_22_", 1:11), paste0("P_D_28_", 1:11))
phage_name_factor = as.factor(phg_order)
phage_feature_binary_l$Var2 = factor(phage_feature_binary_l$Var2, phage_name_factor)

phage_plot = ggplot(phage_feature_binary_l, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_manual(values = color_palette, name = "",  labels = c("contain", "not contain")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none', axis.line.x = element_blank(), axis.line.y = element_blank())
phage_plot
ggsave(paste0(PATH,"fig1_hm_phage_mut.svg"), width=7,height = 7, plot = phage_plot)


phage_feature_binary_add_an_l = melt(phage_only_eqtl_m_allmut_allgen)
phage_feature_binary_add_an_l$value = factor(phage_feature_binary_add_an_l$value, levels = c(0, 1))

row_order_an = rev(sort(unique(phage_feature_binary_add_an_l$Var1)))
row_order_an_factor = as.factor(row_order_an)
phage_feature_binary_add_an_l$Var1 = factor(phage_feature_binary_add_an_l$Var1, row_order_an_factor)

phg_order_an = c("an_p", paste0("P_D_8_", 1:11), paste0("P_D_15_", 1:11), paste0("P_D_22_", 1:11), paste0("P_D_28_", 1:11))
phage_name_an_factor = as.factor(phg_order_an)
phage_feature_binary_add_an_l$Var2 = factor(phage_feature_binary_add_an_l$Var2, phage_name_an_factor)

phage_plot_an = ggplot(phage_feature_binary_add_an_l, aes(Var2, Var1)) + 
  geom_tile(aes(fill = value)) + 
  scale_fill_manual(values = color_palette, name = "",  labels = c("contain", "not contain")) + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"), axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position='none', axis.line.x = element_blank(), axis.line.y = element_blank())

phage_plot_an
ggsave(paste0(PATH,"fig1_hm_phage_mut_anc.svg"), width=7,height = 7, plot = phage_plot_an)
