countsdata<-anti
names(countsdata)<-c("Geneid","D12sense","D4sense",
                     "L12sense","L4sense","D12anti",
                     "D4anti","L12anti","L4anti")
countdata<-countsdata[,2:ncol(countsdata)]
rownames(countdata)<-countsdata[,1]
colnames(countdata)
summary(countdata$D4sense)
summary(countdata$L12sense)
library(ggplot2)
colnames(countdata)
input <- countdata[,c("D12anti","L4anti","L12anti","D4anti")]
print(head(input))
b<-log2(input)
library(reshape2)
library(ggpubr)
b_comparisons <- list(c("L4anti", "D12anti"), c("L12anti", "L4anti"), c("D4anti", "L12anti"),c("D12anti", "D4anti"))
# options(repr.plot.width=4, repr.plot.height=4)
b<-a[,1:4]
data_m<-melt(b)
names(data_m)<-c('sample','log2FPKM')
mytheme2 <- theme(plot.title = element_text(face = "bold.italic",
                                            size = "14", color = "brown"),
                  axis.title = element_text(face = "bold.italic",
                                            size = "10",color = "blue"),
                  axis.text.x = element_text(face = "bold",
                                             size = 8, angle = 45, hjust = 1, vjust = 1),
                  axis.text.y = element_text(face = "bold",size = 8),
                  panel.background = element_rect(fill = "white", 
                                                  color = "black"),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 10,
                                              face = "bold"),
                  panel.grid.minor.x = element_blank())
p <- ggplot(data_m, aes(x=sample, y=log2FPKM),color=sample) + 
        # geom_violin(aes(fill=factor(sample))) + 
        geom_boxplot( aes(fill=factor(sample)),width= .2) +
        # stat_compare_means(comparisons = b_comparisons,label.y = c(10, 11,12,13))+
        # stat_compare_means(label.y = 20, label.x = 1.5)+
        stat_summary(fun.y="mean",geom="point",position = position_dodge(0.8),shape=23,size=1,fill="white")+
        theme(axis.text.x=element_text(angle=50,hjust=0.5, vjust=0.5)) +
        theme(legend.position="none")
a=p+mytheme2
a
# Plot the chart.
boxplot(b, xlab = "D4vsL4_l4",
        ylab = "log2(reads+1)", 
        varwidth = TRUE, 
        col = c("green","yellow"))

h<-countdata[countdata$D4sense<countdata$L4sense, ]
input <- h[,c("D4sense","L4sense")]
print(head(input))
b<-log2(input+1)
# Plot the chart.
boxplot(b, xlab = "D4vsL4_l4UP",
        ylab = "log2(reads+1)", 
        varwidth = TRUE, 
        col = c("gold","red"))


t<-countdata[countdata$D4sense>countdata$L4sense, ]
input <- t[,c("D4sense","L4sense")]
print(head(input))
b<-log2(input+1)
# Plot the chart.
boxplot(b, xlab = "D4vsL4_l4DOWN",
        ylab = "log2(reads+1)", 
        varwidth = TRUE, 
        col = c("gold","red"))
