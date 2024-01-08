library(reshape2)
library(ggplot2)
#setwd("~/Documents/2019coevolution_geno_pheno/code_and_data_paper3/figs_als")

train_res_l = melt(t(all_logistic_res_allf[1:6, ]))
train_res_l$Var1=rep("training",nrow(train_res_l ))
############ FIG S2 A
pdf('supp_fig2_1_bis.pdf', height = 6, width = 4, useDingbats = F)
par(las = 2, mar = c(8, 4, 2, 2))
# boxplot(value ~ Var2, data = train_res_l, names = c("Null", "H:MF", "P:MF", "P+H:MF", "PxH:MF", "Joint:MF"), ylab = "Classification error", main = "Training error on first step logistic regression", col = rgb(0, 0, 0, 0.5), notch = T)
boxplot(value ~ Var2, data = train_res_l,xlab="", xaxt='n', ylab = "Classification error", col = rgb(0, 0, 0, 0.5), notch = T)

tick <- seq_along(c("Null", "H:MF", "P:MF", "P+H:MF", "PxH:MF", "Joint:MF"))
axis(1, at = tick, labels = F)
text(tick, 0.065, c("Null", "H:MF", "P:MF", "P+H:MF", "PxH:MF", "Joint:MF"), srt = 45, xpd = T)

dev.off()


training_lm_res_allf = melt(t(boot_2nd_step_run_mae_res_all[1:6, ]))
training_lm_res_allf$Var1<-rep("training",nrow(training_lm_res_allf))
############ FIG S2 B
pdf('supp_fig2_2_bis.pdf', height = 6, width = 4, useDingbats = F)
par(las = 2, mar = c(8, 4, 2, 2))
# boxplot(value ~ Var2, data = training_lm_res_allf, names = c("Null", "HoMF", "PoMF", "P&HMF", "PCHMF", "CoMF"), ylab = "MAE", main = "Training MAE on second step linear model", col = rgb(0, 0, 0, 0.5), notch = T)
boxplot(value ~ Var2, data = training_lm_res_allf, xlab="", xaxt='n', ylab = "MAE", col = rgb(0, 0, 0, 0.5), notch = T)
tick <- seq_along(c("Null", "H:MF", "P:MF", "P+H:MF", "PxH:MF", "Joint:MF"))
axis(1, at = tick, labels = F)
text(tick, 0.63, c("Null", "H:MF", "P:MF", "P+H:MF", "PxH:MF", "Joint:MF"), srt = 45, xpd = T)

dev.off()

##### Optimized figure 2 - classification error of POA prediction for the validat/training sets




###################
test_res_l$Var1=rep("validation",nrow(test_res_l))
############ FIG 2 A
pdf('fig2_1_als.pdf', height = 6, width = 4, useDingbats = F)
par(las = 2, mar = c(8, 4, 2, 2))
# boxplot(value ~ Var2, data = test_res_l, names = c("Null", "HoMF", "PoMF", "P&HMF", "PCHMF", "CoMF"), ylab = "Classification error", main = "Test error on first step logistic regression", col = rgb(0, 0, 0, 0.5), notch = T)
boxplot(value ~ Var2, data = test_res_l,xlab="", xaxt='n' , ylab = "Classification error", col = rgb(0, 0, 0, 0.5), notch = T)
tick <- seq_along(c("Null", "H:MF", "P:MF", "P+H:MF", "PxH:MF", "Joint:MF"))
axis(1, at = tick, labels = F)
text(tick, 0.065, c("Null", "H:MF", "P:MF", "P+H:MF", "PxH:MF", "Joint:MF"), srt = 45, xpd = T)

dev.off()

#### Combination of fig 2 and 2A
train_test_l<-rbind(train_res_l, test_res_l)
train_test_l$Var2<-gsub("train_","",train_test_l$Var2)
train_test_l$Var2<-gsub("test_","",train_test_l$Var2)
train_test_l$Var2 <- factor(train_test_l$Var2 , levels=c("null", "fof_bo", "fof_po", "fof","sof","cof"))


fig2a=ggplot(train_test_l, aes(x=Var2, y=value, color=Var1, shape=Var1)) + 
  geom_boxplot(        
    # custom boxes
    # color=c("blue","red"),
    # fill=Var1,
    # alpha=0.2,
    # Notch?
    notch=TRUE,
    notchwidth = 0.8,
    
    # custom outliers
    outlier.colour="red",
    outlier.fill="white",
    outlier.size=1,
    lwd=0.5)+
  scale_color_manual(name="set",values = c("grey", "black"))+
  #geom_jitter(size=0.8, alpha=0.5) +
  scale_y_continuous(name = "Classification error")+ theme_classic()+
scale_x_discrete(labels=c("null", "H-only", "P-only", "linear", "non-linear", "mixed"),name ="model")+ theme(legend.position="none", text = element_text(size = 20))
fig2a
ggsave("Fig2A.svg",plot=fig2a,width=5,height=6)


test_lm_res_allf = melt(t(boot_2nd_step_run_mae_res_all[7:12, ]))
test_lm_res_allf$Var1<-rep("validation",nrow(test_lm_res_allf ))
############ FIG 2 B
pdf('fig2_2_als.pdf', height = 6, width = 4, useDingbats = F)
par(las = 2, mar = c(8, 4, 2, 2))
# boxplot(value ~ Var2, data = test_lm_res_allf, names = c("Null", "HoMF", "PoMF", "P&HMF", "PCHMF", "CoMF"), ylab = "MAE", main = "Test set MAE on second step linear model", col = rgb(0, 0, 0, 0.5), notch = T)

boxplot(value ~ Var2, data = test_lm_res_allf, xlab="", xaxt='n'  , ylab = "MAE", col = rgb(0, 0, 0, 0.5), notch = T)
tick <- seq_along(c("Null", "H:MF", "P:MF", "P+H:MF", "PxH:MF", "Joint:MF"))
axis(1, at = tick, labels = F)
text(tick, 0.79, c("Null", "H:MF", "P:MF", "P+H:MF", "PxH:MF", "Joint:MF"), srt = 45, xpd = T)


dev.off()


#### Combination of fig 2B and S2B
train_test_l<-rbind(training_lm_res_allf, test_lm_res_allf)
train_test_l$Var2<-gsub("train_","",train_test_l$Var2)
train_test_l$Var2<-gsub("test_","",train_test_l$Var2)
train_test_l$Var2 <- factor(train_test_l$Var2 , levels=c("null", "fof_bo", "fof_po", "fof","sof","cof"))

fig2b=ggplot(train_test_l, aes(x=Var2, y=value, color=Var1, shape=Var1)) + 
  geom_boxplot(        
    # custom boxes
    # color=c("blue","red"),
    # fill=Var1,
    # alpha=0.2,
    # Notch?
    notch=TRUE,
    notchwidth = 0.8,
    
    # custom outliers
    outlier.colour="red",
    outlier.fill="white",
    outlier.size=1,
    lwd=0.5)+
  scale_color_manual(name="set",values = c("grey", "black"))+
  #geom_jitter(size=0.8, alpha=0.5) +
  scale_y_continuous(name = "Mean absolute error")+ theme_classic()+
  scale_x_discrete(labels=c("null", "H-only", "P-only", "linear", "non-linear", "mixed"),name ="model")+ theme(legend.position="none", text = element_text(size = 20))
fig2b
ggsave("Fig2B.svg",plot=fig2b, width=5,height=6)


##########Statistical differences
library(multcomp)
df<-test_lm_res_allf[,c(2:3)]
df<-test_res_l[,c(2:3)]
res_aov <- aov(df$value ~ Var2,data=df)
summary(res_aov)
#<2e-16
#<2e-16
# Tukey HSD test:
post_test <- glht(res_aov,linfct = mcp(Var2 = "Tukey"))
summary(post_test)
#<1e-04 
#<1e-04 
TukeyHSD(res_aov)
# Groups POA  A=(cof,sof,fof),B=null,C=P,D=H)
#Groups EFF A=(cof,sof),B=null,C=P,D=H,E=fof)

