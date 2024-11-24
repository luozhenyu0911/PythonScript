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
mytheme1 <- theme(plot.title = element_text(face = "bold", hjust = 0.5, vjust = 0.5,
                                            size = "14", color = "black"),
                  axis.title = element_text(face = "bold",
                                            size = "13",color = "black"),
                  axis.text.x = element_text(face = "bold",
                                             size = 13,  hjust = 0.5, vjust = 0.5),
                  axis.text.y = element_text(face = "bold",size = 13),
                  panel.background = element_rect(fill = "white", 
                                                  color = "black"),
                  # axis.line = element_line(size = 1,colour = "red"),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
				  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  legend.text = element_text(size = 8),
                  legend.title = element_text(size = 10,
                                              face = "bold"),
                  panel.grid.minor.x = element_blank())
				  
mytheme3<-theme(axis.title = element_text(face = "bold",
                                            size = "13",color = "black"),
                  axis.text.x = element_text(face = "bold",
                                             size = 13,  hjust = 0.5, vjust = 0.5),
                  axis.text.y = element_text(face = "bold",size = 13)) 
				  
				  
				  
				  
				  
				  
				  
				  
				  
				  

my_bar_theme <- theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 16), 
                      # 因为 x 轴标签要旋转 90°，所以这里用来旋转
                      axis.text.y = element_text(size = 16),
                      axis.title.y = element_text(size = 18,
                                                  face = "bold",),
                      legend.position = "none") +
  theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5)) # 基因名要居中，这里用来居中。
