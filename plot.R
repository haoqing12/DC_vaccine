library(tidyverse)
library(data.table)
library(swimplot)

p<-df %>%
  ggplot(aes(y=Patients_ID, group=Patients_ID)) +
  theme_bw() +
  geom_line(aes(x=value, col = variable),,size=2)+
  scale_color_manual(values = c( "#74aaff","#483d8b","#F4CD57","#CCCCCC","#b24745",'red'))+
  ylab(" ")+
  xlab("")+
  scale_x_continuous(expand=c(0,0))


p<-ggplot(df1, aes(x=variable, y= value1, fill = phase))+ 
    geom_jitter(aes(fill = phase, color = phase),
        position = position_jitterdodge(dodge.width = 0.3,jitter.width = 0.5),
        shape = 21,size = 1)+ 
  geom_boxplot(aes(group = variable,color = phase), 
        width = 0.2,fill = "white",
        outlier.color = NULL,outlier.alpha =0,
        outlier.fill = NULL,outlier.shape = NULL)+
  scale_color_manual(values = c("#EB3323",  "#5BBCD6"))+
  scale_fill_manual(values = c("#EB3323",  "#5BBCD6"))+
  ylab("IFNγ+ spots")+
  xlab("")+
  facet_grid(cols = vars(phase))+
  theme_bw()+
    scale_y_continuous(expand = c(0, 0)) +
    geom_hline(aes(yintercept=sqrt(20)),color = 'grey', linetype="dashed")+
    geom_hline(aes(yintercept=sqrt(200)),color = 'grey', linetype="dashed")


p<-ggplot(df3, aes(x=variable, y= 100*value))+ 
  geom_boxplot(aes(group = variable,color = phase), 
        width = 0.3,fill = "white",
        outlier.color = NULL,outlier.alpha =0,
        outlier.fill = NULL,outlier.shape = NULL)+
    geom_jitter(aes(fill = phase, color = phase), 
    position = position_jitterdodge(dodge.width = 0.3,jitter.width = 0.5),
    shape = 21,size = 3, color = "white")+ 
    scale_fill_manual(values =  c("#EB3323",  "#5BBCD6"))+
    scale_color_manual(values = c("#EB3323",  "#5BBCD6"))+
  ylab("Positive rate %")+
  xlab("")+
  facet_grid(cols = vars(phase))+
  theme_bw()+
  stat_compare_means(aes(label = paste0("p = ",..p.format..), group = "phase"), method = "wilcox.test")+
    scale_y_continuous(expand = c(0, 0))+
    guides(fill = "none") 


p<-ggplot(df4 , aes(x=ID, y= post, fill = Trunk_driver))+ 
#   geom_point(aes(color = Trunk_driver, group = ID), size=1) +
 geom_jitter(
        aes(fill = Trunk_driver),
        # position = position_dodge(0.3),
        position = position_jitterdodge(
        # 类似于position_dodge
        dodge.width = 0.3,
        #类似于bar的宽度
        jitter.width = 0.2),
        shape = 21,
        color = "black",
        size = 2.5)+ 
   ylab("IFNγ+ spots of post-vaccination")+
  xlab("")+
  theme_bw()+
  scale_fill_manual(values = c("white",  "red", "#F7D7B3"))+
    scale_y_continuous(expand = c(0, 0)) 