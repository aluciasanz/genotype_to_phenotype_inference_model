options(stringsAsFactors = F)
library(pracma)
PATH='./figures_04252022/'

list_coef_phage_plot <- function(phage_coef_list, min_v, max_v){
  plotobj = ggplot(phage_coef_list, aes(x=x, y=y)) + 
    geom_tile(aes(x=x, y=y, fill = coef), col = "white") + 
    geom_segment(aes(x = 0.5, xend = 178.5, y = 0.5, yend = 0.5, size = 0.1)) + 
    geom_segment(aes(x = 0.5, xend = 178.5, y = 1.5, yend = 1.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 0.5, y = 0.5, yend = 1.5, size = 0.2)) + 
    geom_segment(aes(x = 178.5, xend = 178.5, y = 0.5, yend = 1.5, size = 0.2)) + 
    scale_size_identity() + 
    xlab("Phage Mutations") + 
    ylab("Weights") +
    scale_fill_gradientn(colours = c("darkred", "red", "white","green", "darkgreen"), values = scales::rescale(c(min_v, max_v)), guide = "colorbar", limits=c(min_v, max_v), breaks=linspace(min_v, max_v, 5)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), panel.border = element_rect(fill = NULL, color = 'white'))
  return(plotobj)
}

list_coef_bac_plot <- function(bac_coef_list, min_v, max_v){
  plotobj = ggplot(bac_coef_list, aes(x=x, y=y)) + 
    geom_tile(aes(x=x, y=y, fill = coef), col = "black", size = 0.2) + 
    xlab("Bac Mutations") + 
    ylab("Weights") +
    scale_fill_gradientn(colours = c("darkred", "red", "white","green", "darkgreen"), values = scales::rescale(c(min_v, max_v)), guide = "colorbar", limits=c(min_v, max_v), breaks=linspace(min_v, max_v, 5)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), panel.border = element_rect(fill = NULL, color = 'white'))
}

matrix_coef_plot <- function(mut_mat_l, min_v, max_v){
  ggplot(mut_mat_l, aes(x=phage_mut, y=bac_mut)) + 
    geom_tile(aes(x=phage_mut, y=bac_mut, fill = coef), col = "white") + 
    geom_segment(aes(x = 0.5, xend = 178.5, y = 0.5, yend = 0.5, size = 0.2)) + 
    geom_segment(aes(x = 0.5, xend = 178.5, y = 1.5, yend = 1.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 2.5, yend = 2.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 3.5, yend = 3.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 4.5, yend = 4.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 5.5, yend = 5.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 6.5, yend = 6.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 7.5, yend = 7.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 8.5, yend = 8.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 9.5, yend = 9.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 10.5, yend = 10.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 11.5, yend = 11.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 12.5, yend = 12.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 13.5, yend = 13.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 14.5, yend = 14.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 15.5, yend = 15.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 16.5, yend = 16.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 17.5, yend = 17.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 18.5, yend = 18.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 19.5, yend = 19.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 0.5, y = 0.5, yend = 19.5, size = 0.2)) + 
    geom_segment(aes(x = 178.5, xend = 178.5, y = 0.5, yend = 19.5, size = 0.2)) + 
    scale_size_identity() +
    xlab("Phage Mutations") + 
    ylab("Host Mutations") +
    scale_x_continuous(limits = c(0, 181)) + 
    scale_y_reverse(limits = c(20, 0)) +
    scale_fill_gradientn(colours = c("darkred", "red","white","green", "darkgreen"), values = scales::rescale(c(min_v, max_v)), guide = "colorbar", limits=c(min_v, max_v), breaks=linspace(min_v, max_v, 5)) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank(), axis.line.x = element_blank(), axis.line.y = element_blank(), panel.border = element_rect(fill = NULL, color = 'white'))
}


matrix_coef_plot_unused <- function(mut_mat_l, min_v, max_v){
  ggplot(mut_mat_l, aes(x=phage_mut, y=bac_mut)) + 
    geom_tile(aes(x=phage_mut, y=bac_mut, fill = coef), col = "black") + 
    geom_segment(aes(x = 0.5, xend = 178.5, y = 0.5, yend = 0.5, size = 0.2)) + 
    geom_segment(aes(x = 0.5, xend = 178.5, y = 1.5, yend = 1.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 2.5, yend = 2.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 3.5, yend = 3.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 4.5, yend = 4.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 5.5, yend = 5.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 6.5, yend = 6.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 7.5, yend = 7.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 8.5, yend = 8.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 9.5, yend = 9.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 10.5, yend = 10.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 11.5, yend = 11.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 12.5, yend = 12.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 13.5, yend = 13.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 14.5, yend = 14.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 15.5, yend = 15.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 16.5, yend = 16.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 17.5, yend = 17.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 18.5, yend = 18.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 178.5, y = 19.5, yend = 19.5, size = 0.2)) +
    geom_segment(aes(x = 0.5, xend = 0.5, y = 0.5, yend = 19.5, size = 0.2)) + 
    geom_segment(aes(x = 178.5, xend = 178.5, y = 0.5, yend = 19.5, size = 0.2)) + 
    scale_size_identity() +
    xlab("Phage Mutations") + 
    ylab("Host Mutations") +
    scale_x_continuous(limits = c(0, 180)) + 
    scale_y_reverse(limits = c(19, 0)) +
    scale_fill_gradientn(colours = c("darkred", "red","white","green", "darkgreen"), values = scales::rescale(c(min_v, max_v)), guide = "colorbar", limits=c(min_v, max_v), breaks=linspace(min_v, max_v, 5))
}

# 1st step fof bac only
final_1st_fof_bo_min
final_1st_fof_bo_max
final_1st_fof_bo_bac = list_coef_bac_plot(first_step_fof_bo_bac_beta_l, -5.5, 5.5)
final_1st_fof_bo_bac
ggsave(paste0(PATH,"supp_fig6_fof_bo_bac_coef.pdf"), plot = final_1st_fof_bo_bac)

# 1st step fof phage only
final_1st_fof_po_min
final_1st_fof_po_max
final_1st_fof_po_phage = list_coef_phage_plot(first_step_fof_po_phage_beta_l, -2.5, 2.5)
final_1st_fof_po_phage 
ggsave(paste0(PATH,"supp_fig6_fof_po_phage_coef.pdf"), plot = final_1st_fof_po_phage)

# 1st step fof
final_1st_fof_min
final_1st_fof_max
final_1st_fof_phage = list_coef_phage_plot(first_step_fof_phage_beta_l, -8.5, 8.5)
final_1st_fof_bac = list_coef_bac_plot(first_step_fof_bac_beta_l, -8.5, 8.5)

ggsave(paste0("fig3_fof_phage_coef.svg"), plot = final_1st_fof_phage)
ggsave(paste0(PATH,"fig3_fof_bac_coef.pdf"), plot = final_1st_fof_bac)

# 1st step sof, matrix order is row from up to down 1-18, left to right 1-176
final_1st_sof_min
final_1st_sof_max
final_1st_sof_m = matrix_coef_plot(first_step_sof_beta_matrix_l, -35.5, 35.5)

ggsave(paste0("fig3_sof_m_coef.pdf"), plot = final_1st_sof_m)

# 1st step cof
final_1st_cof_min
final_1st_cof_max
final_1st_cof_phage = list_coef_phage_plot(first_step_cof_phage_beta_l, -18.6, 18.6)
final_1st_cof_bac = list_coef_bac_plot(first_step_cof_bac_beta_l, -9, 9)
final_1st_cof_m = matrix_coef_plot(first_step_cof_beta_matrix_l, -35.5, 35.5)

ggsave(paste0(PATH,"fig3_cof_phage_coef.pdf"), plot = final_1st_cof_phage)
ggsave(paste0(PATH,"fig3_cof_bac_coef.pdf"), plot = final_1st_cof_bac)
ggsave(paste0("fig3_cof_m_coef.pdf"), plot = final_1st_cof_m)


# 2nd step fof bac only
final_2nd_fof_bo_min
final_2nd_fof_bo_max
final_2nd_fof_bo_bac = list_coef_bac_plot(second_step_fof_bo_bac_beta_l, -2.1, 2.1)

ggsave(paste0(PATH,"supp_fig7_fof_bo_bac_coef.svg"), plot = final_2nd_fof_bo_bac)

# 2nd step fof phage only
final_2nd_fof_po_min
final_2nd_fof_po_max
final_2nd_fof_po_phage = list_coef_phage_plot(second_step_fof_po_phage_beta_l, -3.6, 3.6)

ggsave(paste0(PATH,"supp_fig7_fof_po_phage_coef.svg"), plot = final_2nd_fof_po_phage)

# 2nd step fof
final_2nd_fof_min
final_2nd_fof_max
final_2nd_fof_phage = list_coef_phage_plot(second_step_fof_phage_beta_l, -3, 3)
final_2nd_fof_bac = list_coef_bac_plot(second_step_fof_bac_beta_l, -3, 3)

ggsave(paste0(PATH,"fig4_fof_phage_coef.svg"), plot = final_2nd_fof_phage)
ggsave(paste0(PATH,"fig4_fof_bac_coef.svg"), plot = final_2nd_fof_bac)

# 2nd step sof
final_2nd_sof_min
final_2nd_sof_max
final_2nd_sof_m = matrix_coef_plot(second_step_sof_beta_matrix_l, -3, 3)

ggsave(paste0("fig4_sof_m_coef.pdf"), plot = final_2nd_sof_m)

# 2nd step cof
final_2nd_cof_min
final_2nd_cof_max
final_2nd_cof_phage = list_coef_phage_plot(second_step_cof_phage_beta_l, -3, 3)
final_2nd_cof_bac = list_coef_bac_plot(second_step_cof_bac_beta_l, -3, 3)
final_2nd_cof_m = matrix_coef_plot(second_step_cof_beta_matrix_l, -3, 3)

ggsave(paste0("fig4_cof_phage_coef.pdf"), plot = final_2nd_cof_phage)
ggsave(paste0("fig4_cof_bac_coef.pdf"), plot = final_2nd_cof_bac)
ggsave(paste0("fig4_cof_m_coef.pdf"), plot = final_2nd_cof_m)