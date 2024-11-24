
data<-shijiain
library(ggplot2)
colnames(data)
colnames(shijiain)
ggplot(data, aes(x = "regulated", y ="DEG number" , fill = factor(Type) ,
                            order = as.numeric(factor(data$Type,
                                                      levels=c("+","-","o")))   )) + 
      geom_bar(stat = "identity", width = 1, position = "stack" ) +  
      facet_grid(. ~ Name,scales = "free_x") +
       theme(axis.text.x =element_blank())+mytheme1
p
p1<-p1+mytheme1

channel2<-shijian
p<-ggplot(data = channel2, aes(x = Left, y = Amount, label =Amount ,fill = Type ),
                            order = as.numeric(factor(channel2$Type,
                                   levels=c("*","-","+","*","-","+")))   )+
  geom_bar(stat = "identity", width = 1, position = "stack" ) +  
  facet_grid(cols = vars(Name),scales = "free_x") +
  theme(axis.text.x =element_blank())+mytheme1+
  geom_text(size = 3, position = position_stack(vjust = 0.5))
p  
p1<-p+labs(x = "",
           y = "Sense DEG Amount")  
p1
