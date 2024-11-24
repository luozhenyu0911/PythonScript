rm(list=ls())
library(data.table)
library(ggplot2)
library(ggsignif)
colnames(geneid)
Group1<-X1geneidlog
Group2<-geneidlog
# 整合数据到一个数据框：
b<-rbind(Group1,Group2)
b<-melt(b,id.vars = c("group"))
b$group<-as.factor(b$group)
c<-copy(b)
setDF(c)
c1<-tapply(c[c$group==1,"value"],c[c$group==1,"variable"],mean)
c2<-tapply(c[c$group==2,"value"],c[c$group==2,"variable"],mean)
c3<-rbind(data.frame(variable=names(c1),value=c1,group=1),data.frame(variable=names(c2),value=c2,group=2))

c3$group<-as.factor(c3$group)
# 分别计算两组均值，用来画折线图：
c3$variable2<-NA
c3[c3$group==1&c3$variable=="D4","variable2"]<-0.795
c3[c3$group==1&c3$variable=="D12","variable2"]<-1.795
c3[c3$group==1&c3$variable=="L4","variable2"]<-2.795
c3[c3$group==1&c3$variable=="L12","variable2"]<-3.795


c3[c3$group==2&c3$variable=="D4","variable2"]<-1.185
c3[c3$group==2&c3$variable=="D12","variable2"]<-2.185
c3[c3$group==2&c3$variable=="L4","variable2"]<-3.185
c3[c3$group==2&c3$variable=="L12","variable2"]<-4.185


p1<-ggplot(b)+
  geom_boxplot(aes(x=variable,y=value,fill=group),width=0.6,position = position_dodge(0.8),outlier.size = 0,outlier.color = "white")+
  scale_fill_manual(values = c("red", "Blue"),breaks=c("1","2"),labels=c("Group 1","Group 2"))+
  geom_point(data=c3,aes(x=variable2,y=value,color=group),shape=15,size=2,fill="white")+
  geom_line(data=c3,aes(x=variable2,y=value,color=group),size=1.5,linetype = "dotted")+
  # geom_smooth(data=c3,aes(x=variable2,y=value,color=group),size=1,linetype = "dashed")+
  # stat_summary(fun.y = mean, geom = "errorbar", aes(x=variable,y=value,ymax = ..y.., ymin = ..y..,color=group),width = .75, linetype = "dashed")
  xlab("sample")+
  ylab("Standardized_expression")+
  theme_bw()+
  theme(
    legend.position = "top",
    legend.background=element_blank(),
    legend.key = element_blank(),
    legend.margin=margin(0,0,0,0,"mm"),
    axis.text.x=element_text(size=rel(1.1),face="bold"),
    axis.line.x = element_line(size = 0.5, colour = "black"),
    axis.line.y = element_line(size = 0.5, colour = "black"),
    legend.text=element_text(size=rel(1.1)),
    legend.title=element_blank(),
    panel.border = element_blank(),
    panel.grid = element_blank()
  )+
  guides(color=FALSE)
p2<-p1+mytheme2
p2
