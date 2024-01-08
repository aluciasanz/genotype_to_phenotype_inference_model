options(stringsAsFactors = F)
library(ridge)
library(Matrix)
library(gridExtra)
#for linux you may need apt install cmake AND package nloptr
library(ggpubr)
#setwd("~/Documents/2019coevolution_geno_pheno/code_and_data_paper3/figs_als")
#dist on eop sup figure
PATH='./figures_04252022/'
sptest = shapiro.test(log(label_allgen[label_allgen > 0]))
log_pos_label = log(label_allgen[label_allgen > 0])
#Final figure for EOF
histeop<-ggplot(data.frame(x=data[data>0]),aes(x=x))+
  geom_histogram(color="black",
                 fill="white",
                 lwd=0.5, 
                 position = "identity",bins=101,breaks=seq(0,5,by=0.1))+
  xlab("EOP>0")+
  theme_classic()+ theme( text = element_text(size = 20))
ggsave("FigS1_inset.svg",plot=histeop,width=5,height = 5)
#
pdf(paste0(PATH,"supp_fig13_eop_distribution.pdf"), height = 10, width = 10, useDingbats = F)
par(mfrow = c(2,2))
hist(log_pos_label, breaks = 40, col = 'gray', main = "Distribution of log positive EOP values", xlab = "log(EOP)")
dev.off()
pdf(paste0(PATH,"supp_fig1_1.pdf"), height = 10, width = 10, useDingbats = F)
hist(label_allgen, breaks = 100, col = 'gray', main = "Distribution of all EOP values", xlab = "EOP value")
dev.off()
pdf(paste0(PATH,"supp_fig1_2.pdf"), height = 10, width = 10, useDingbats = F)
hist(label_allgen[label_allgen > 0], breaks = 40, col = 'gray', main = "Distribution of positive EOP values only", xlab = "EOP value")
dev.off()

pdf(paste0(PATH,"supp_fig3_qqplot_2.pdf"), height = 5, width = 5, useDingbats = F)
ggqqplot(log_pos_label)
dev.off()
pdf(paste0(PATH,"supp_fig3_qqplot_1.pdf"), height = 5, width = 5, useDingbats = F)
hist(log_pos_label, breaks = 40, main='',col = 'gray', xlab = "log(EOP)")
dev.off()