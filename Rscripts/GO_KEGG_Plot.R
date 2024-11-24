#GO多组的热图见 tidyr::spread.R。
rm(list = ls())
options(stringsAsFactors = F)

library(XLConnect)
library(ggplot2)
#library(cowplot)
library(export)

#---------------------------单独上下调的GO bubble plot
#for (i in c(1:9)){

# x轴为GeneRatio,适合单个DEGs的GO富集的情况；x轴为Group分组,适合多组DEGs的GO富集的情况。
egoup<- readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',4)
nrow(egoup)

#按照GeneRatio列进行排序
order<-order(egoup$GeneRatio)
egoup$GOterms<- factor(egoup$GOterms, 
                       levels = unique(egoup$GOterms[order]))

cc = colorRampPalette(c('red','white'))
col=cc(100)[5:50]   #去除极颜色
col=rev(col)

gobubbleup <- ggplot(egoup,aes(GeneRatio,GOterms)) +
  theme_bw(base_size = 18) +#ggtitle('CK-96h/24h-down')+
  geom_point(aes(size = Counts,color=-log10(FDR)),alpha=1)+  # color=FDR or padj
  #scale_colour_gradient(low="red",high="blue")+
  scale_color_gradientn(colours=col)+
  scale_size_continuous(range = c(5,10))+
  labs(#size="Gene counts",  #colour=expression(-log10("P Value")),,title=main.title
       x="Gene Ratio",y="")+
  theme(title =element_text(size = 18) ,
        axis.text.x  = element_text(size = 14,colour = 'black'),
        axis.text.y  = element_text(size = 18,colour = 'black'),
        axis.title.x = element_text(size = 18,face = 'bold'), 
        axis.title.y = element_text(size = 18),
        #panel.background = element_rect(fill='transparent',color ="black"),
        #panel.border =element_rect(fill='transparent')
        )+
  facet_grid(Ontology~.,scales= "free",space = "free")
gobubbleup

graph2ppt(file='4clusters-GO-BP.pptx')


# x轴为Group分组,适合多组DEGs的GO富集的情况。x轴为GeneRatio,适合单个DEGs的GO富集的情况。
egogroups<- readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',3)
nrow(egogroups)

egogroups<-egogroups[order(egogroups$FDR,decreasing = T),]
egogroups<-egogroups[order(egogroups$GeneRatio,decreasing = F),]# or Counts
egogroups$GOterms<-factor(egogroups$GOterms,
                      levels=unique(as.character(egogroups$GOterms)))

cc = colorRampPalette(c('red','white'))
col=cc(100)[1:60]   #去除极颜色
col=rev(col)

gobubble <- ggplot(egogroups,aes(Group,GOterms))+
  theme_bw(base_size = 18)+#ggtitle('CK-96h/24h-down')+
  geom_point(aes(size = Counts,color=-log10(FDR)),
             alpha=1)+  # color=FDR or padj
  scale_color_gradientn(colours=col)+
  #scale_color_viridis_c(begin = 0.3, end = 1) +  ## 调整颜色的区间,begin越大，整体颜色越明艳
  scale_size_continuous(range = c(5,10))+
  labs(x="",y="")+
  theme(legend.position = c(-2,0.5),
        legend.title = element_text(face = 'bold'),
        strip.text = element_text(face = 'bold'),
        title =element_text(size = 18) ,
        axis.text.x  = element_text(size = 18,#face = 'bold',
                                    colour = 'black'),
        axis.text.y  = element_text(size = 18,#face = 'bold',
                                    colour = 'black'),
        axis.title.x = element_text(size = 0), 
        axis.title.y = element_text(size = 0),
        #panel.background = element_rect(fill='transparent',color ="black")
        )  +  facet_grid(Ontology~.,scales= "free",
                         space = "free")

gobubble

graph2ppt(file='GO-BP.dotplot.pptx')
dev.off()
#}

#合并两图
library(ggpubr)
p <- ggarrange(gobubbleup+rremove("xlab"), gobubbledown+rremove("ylab")+rremove("xlab"),
               common.legend = TRUE,nrow = 1, legend = "right")
figure<- annotate_figure(p,top = text_grob("48h-LN_CK", color = "black", face = "bold", size = 14),
                         bottom = text_grob("GeneRatio",hjust = 0.5),
                         #fig.lab = "GeneRatio",fig.lab.pos  = 'bottom',fig.lab.size =12
                         )
figure

#plot_grid(gobubbleup,gobubbledown,labels = c('Up','Down'))

#--------------------------------------------------------------barplot, 对MF,BP,CC按照颜色展示
rm(list = ls())
options(stringsAsFactors = F)

library(ggpubr)
egoup1<- readWorksheetFromFile('e:/projects/rice/strand_specific/go/GO.xlsx',1)[1:5,]
ggbarplot(egoup1, x = "GOterms", y = "log",
          fill = "Ontology",               # change fill color by Ontology
          color = "white",            # Set bar border colors to white
          palette = "jco",            # 可自定义，jco journal color palett. see ?ggpar
          sort.val = "asc",           # Sort the value in dscending order
          sort.by.groups = TRUE,      # 按分组内进行排序
          #x.text.angle = 90  
          ) + # Rotate vertically x axis texts
  coord_flip()  +labs(
    #size="Gene counts",  #colour=expression(-log10("P Value")),,title=main.title
    x="",y="-log10(FDR)")+ ggtitle('24h-LN/CK-up')+ 
  theme_classic(base_size = 18)


#----------------------------------------------------------------------Bar plot,上下调一起展示
#GO的barplot：x=Description, y = padj
# 统计
up<-readWorksheetFromFile('GO.xlsx',9)[1:30,]
down<-readWorksheetFromFile('GO.xlsx',9)[31:38,]
#合并
go_class <- rbind(up, down)
go_class
go_class$significant <- rep(c("up","down"),times=c(nrow(up),nrow(down)))
#画图
class_go <- ggplot(go_class,aes(x=GOterms, y = -log10(FDR))) +ggtitle('LN-96h/48h')+
  geom_bar(aes(fill=significant),stat="identity",position = "dodge") + 
  scale_fill_manual(values = c(up="#ca0020",down="#2b83ba")) +
  theme(title = element_text(size=16),
        axis.text = element_text(size = 12,colour = 'black'),
        axis.title.x = element_text(size = 14), 
        axis.title.y = element_text(size = 14))+
  theme_bw()  + coord_flip()+ labs(x="-log10(FDR)",y="")  #coord_flip() 坐标轴转置
class_go



##------------------------------------------------------------使用ggplot实现GOplot作图
#适用于 单个 DEGs的GO富集的情况。
library(RColorBrewer)
library(ggrepel)

cc<-rev(brewer.pal(n=7,name = 'RdYlBu'))#配色

ego<-readWorksheetFromFile('ssRNAseq/GO_DEGs.xlsx',2)
#分面
g<-ggplot(ego,aes(x=GeneRatio,y=-log10(FDR)))+
  geom_point(aes(size=Counts,color=Ontology),alpha=1)+
  scale_size(range=c(15,30))+
  scale_color_brewer(palette = 'RdYlBu')+# 配色可选：Greens、Accent、Set1 或数字
  #自行调色,?RColorBrewer::brewer.pal()  或 ?scales::brewer_pal for more details
  theme_bw()+#theme_bw(legend.position = c("none"))
  theme(legend.position = c("none") )+#legend.position = c("none")
  geom_label_repel(data = ego[-log10(ego$FDR)> 0,],#阈值调节
                   aes(label = GOterms,fill=Ontology),size = 3,segment.color = "black",
                   show.legend = FALSE ) +
  facet_grid(.~Group) #以Ontology分面
g
ggsave(g,filename = 'go2.png')



#==============================================================================================GO  plot
#GOplot使用了zscore概念，其并不是指Z-score标准化，而是用于表示参与某个GO Term下基因的上下调情况。

#--------------------------------------------------------------------------------内置数据
library(GOplot)
data(EC)
head(EC$david)#Category    GOID     GOTerm    Genes     adj_pval
head(EC$genelist)# geneID    logFC   AveExpr   t    P.Value      adj.P.Val     B
circ <- circle_dat(EC$david, EC$genelist)
GOBubble(circ, labels = 4) 
#labels: Sets a threshold for the displayed labels. The threshold refers to the -log(adjusted p-value) (default=5)
GOBubble(circ, title = 'Bubble plot', colour = c('orange', 'darkred', 'gold'), 
         display = 'multiple', labels = 3)
#--------------------------------------------------------------------------------------------

#------------------------------------------------------------详细用法说明
#Goplot有 circle_dat函数可以用于整合外来数据，
#但是要求数据的格式（尤其指列名）规范化，可看函数说明或者EC示例文件，模仿着修改。
#文件格式依次包括：category，即三个本体BP,MF,CC   ID:GoID  term:GOterm描述   count   logFC   adj_pval   zscore
library(GOplot)
deg <- read.table(file = "c:/Users/ZQQ/Desktop/deg.txt", #包括两列：差异基因id、logFC,用于计算zscore.
                  sep = "\t", header = T, stringsAsFactors = F)
enrich <- read.table(file = "c:/Users/ZQQ/Desktop/go.txt", #包括category  GoID  term  genes   adj_pval  
                     sep = "\t", header = T, stringsAsFactors = F)
circ <- circle_dat(enrich, deg) #自行整合为所需格式

#barplot图，纵坐标是adj.pvalue，横纵标是GO id，通过柱子颜色来展示zscore
GOBar(circ, display = 'multiple')

#bubble图，以adj.pvalue为纵坐标，zscore为横纵标，
#对于满足一定阈值（比如log10(adj.pvalue) > 3）的显示GO id。
GOBubble(circ, display = "single", labels = 3, table.legend = F)
GOBubble(circ, display = "multiple", labels = 3, table.legend = F)

#GOplot为了过滤一些比较相似的GO ID（认为是冗余的），其根据不同GO之间参与基因的相似程度进行的过滤，
#这个可以自己写代码预先过滤下，也可以直接用其写好的函数reduce_overlap
reduced_circ <- reduce_overlap(circ, overlap = 0.75)
GOBubble(reduced_circ, labels = 3, table.legend = F,display = 'multiple')

#Circular图，内圈展示zscore，外圈展示GO id，外圈里面则是参与基因的上下调分布情况
GOCircle(circ, nsub = 10, label.size = 5, rad1 = 3, rad2 = 4, table.legend = F)

#用圈图展示gene和go term之间的关系，这需要先用chord_dat函数构建一个对象，
#其需要circ，deg以及制定的go id等数据；limit参数则是限制基因/GO的数目.
process <- tail(enrich$ID, 10)
chord <- chord_dat(data = circ, genes = deg, process = process)
GOChord(chord, space = 0.02, gene.order = 'logFC', gene.space = 0.25, gene.size = 2.5, 
        limit = c(0,2))

#热图，其以热图展示GO term和gene id，每行是一个GO id，每列是一个gene；
#如果其图例是count的话，其每格的颜色则是该gene参与图中所有go term总数；
#如果是logFC的话，则比较好理解，是该基因的logFC的值.
GOHeat(chord, nlfc = 1, fill.col = c('red', 'yellow', 'green'))

#圈图，其内圈是根据基因logFC做的聚类，外圈则是GO term的堆积图.
GOCluster(circ, process, clust.by = 'logFC', term.width = 2)
