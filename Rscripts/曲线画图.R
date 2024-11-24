library(ggplot2)
library(reshape2)
library(dplyr)
library(plyr)
new<-atgdegaverage
head(new)
a<-new[,2:5]
rownames(a)<-new$geneid
na<-data.frame(t(scale(t(a))))
na$geneid<-new$geneid
data <- melt(na,id.vars = c('geneid')) #保留id.vars，其余整合。#Cluster  gene  variable  value
data$variable<- factor(data$time,levels= c('D4','D12','L4','L12'),ordered = TRUE)
table(data$variable)
table(data$geneid)
p2<-ggplot(data,aes(x=variable, y=value, color=geneid,group =geneid)) + theme_bw(base_size = 14) + 
  geom_line(size=0.8,alpha=1) +   xlab("")+ ylab("FPKM z-score")+
  # geom_hline(yintercept =0,linetype=2) + #ggtitle("Cluster") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text = element_text(size=12, face = "bold"),
        strip.text = element_text(size = 12, face = "bold"),
  )
p2
ggsave(p1,filename = 'CK-expr-trend.png')
ggsave(p2,filename = 'SiATG')
dev.off()




library(ggplot2)
 library(reshape2)
 n<-a[,1:4]
 # nrow(dfidfm)
# n<-c3[,2:5]
# rownames(n)<-c3[,1]
data<-t(n)
sample<-c("D4","D12","L4","L12")
data<-cbind(sample,data)
data<-data.frame(data,stringsAsFactors = F)
data<-as.data.frame(lapply(data,as.numeric))
data$sample<-c("T1Dark4h","T2Dark12h","T3Light4h","T4Light12h")
# rownames(data)<-rownames(data1)
dfidfm <- melt(data, id.vars="sample")
colnames(dfidfm)
library(dplyr)
h12.coup1 <- dfidfm %>% 
  # filter(Cat == "Cat.1" & group2 == "co_up") %>% 
  group_by(sample) %>% 
  summarise(value=mean(value, na.rm = T))
p<-ggplot(dfidfm, aes(x=sample, y=value),color="blue") + 
      geom_line(aes(group = variable),color="BlueViolet",show.legend = F, alpha=0.15)+
  geom_line(data = dfidfm, aes(group=value), alpha=1,size=1.5)+
  geom_line(data = h12.coup1,aes(group =1), alpha=1,size=1.5)+
  xlab("Treat")+
  ylab("Standardized FPKM")+
  ggtitle("Cluster3")  
p1<-p+mytheme1
p1
   # theme(axis.text.y=element_blank())


targt<-n[which(rownames(n)%in%c("Seita.1G236100.v2.2","Seita.6G055700.v2.2")),]

