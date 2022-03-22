library(ggplot2)
library(ChIPseeker)
library(dplyr)
library(ggrepel)
library(magrittr)
require(ChIPseeker)
require(GenomicFeatures)

#########peak annotation

args=commandArgs(T)
Db <- loadDb(args[1])
x1 = annotatePeak(args[2], tssRegion=c(-2000,0),level="gene",overlap="all",TxDb=Db,assignGenomicAnnotation = TRUE,
                  genomicAnnotationPriority = c("Promoter","5UTR", "3UTR", "Exon", "Intron","Downstream", "Intergenic"))
num <- nrow(read.table(args[2]))
write.table(x1,args[3],sep="\t",row.names = F)

#下游2kb以外的替换为Intergenic
as.data.frame(x1)$annotation %>%
 gsub("^(Intron).*",'\\1',.) %>%
 gsub("^(Exon).*",'\\1',.) %>%
 gsub("Downstream \\(2\\-3kb\\)","Distal Intergenic",.) %>%
 gsub("Downstream \\(1\\-2kb\\)","Downstream \\(<=2kb\\)",.) %>%
 gsub("Downstream \\(<1kb\\)","Downstream \\(<=2kb\\)",.) %>%
 gsub("Promoter \\(1\\-2kb\\)","Promoter \\(<=2kb\\)",.) %>%
 gsub("Promoter \\(<=1kb\\)","Promoter \\(<=2kb\\)",.) %>%
 table(.) %>% data.frame(.)-> peak_sum

write.table(peak_sum,args[4],sep="\t",quote = F,row.names = F,col.names = F)
v<-peak_sum$Freq
g<-peak_sum$.
l<-peak_sum[order(peak_sum$.,decreasing = TRUE),]$.

#########主题
theme_z <- function(..., bg='white'){
  require(grid)
  theme_classic(...) +
    theme(rect=element_rect(fill=bg),
          panel.grid=element_blank(),
          legend.title=element_blank(),
          plot.title=element_text(hjust=0.5,size=15),
          legend.key=element_rect(fill='transparent', color='transparent'),
          legend.text = element_text(color = "black",size=18),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),
          axis.line.y=element_blank())
}
#########绘饼图
df1 <- data.frame(value = v,
                  Group = g) %>%
  mutate(Group = factor(Group, levels = l),
         cumulative = cumsum(value),
         midpoint = cumulative - value / 2,
         label = paste0(round(value / sum(value) * 100,1), "%"))

title<-gsub("_"," ",args[5])
p1<-ggplot(df1, aes(x = 1, weight = value, fill = Group)) +
  geom_bar(width = 1, position = "stack") +
  coord_polar(theta = "y")+theme_z()+
  geom_text_repel(aes(x = 1.7, y = midpoint, label = label),
                  size =8,fontface="bold",color = "black")+
  labs(title = paste(title,num,sep='\n'))+
  scale_fill_manual(values=c("Promoter (<=2kb)"="#2980b9",
                             "5' UTR"="#9b59b6",
                             "3' UTR"="#e74c3c",
                             "Exon"="#16a085",
                             "Intron"="#f1c40f",
                             "Downstream (<=2kb)"="#e67e22",
                             "Distal Intergenic"="#34495e"),
                    limits=c("Promoter (<=2kb)",
                             "5' UTR",
                             "3' UTR",
                             "Exon",
                             "Intron",
                             "Downstream (<=2kb)",
                             "Distal Intergenic"))+
  theme(plot.title=element_text(size=25),
        legend.position = "right",
        legend.justification="center",
        #legend.spacing.y = unit(2.0, 'cm'),
        legend.key = element_rect(size = 8),
        legend.key.height = unit(1, "cm"),
        legend.key.width = unit(1, "cm"))
ggsave(args[6],p1,width = 7,height = 7,units = "in",device = "pdf")

#values=c("Promoter"="#00a4c5","5' UTR"="#eb6e95","3' UTR"="#7dcdf4","Exon"="#ea5532","Intron"="#f0e90f","Downstream"="#f8b500","Distal Intergenic"="#009a53")
