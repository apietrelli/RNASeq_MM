### R script for STAR alignment data Plot
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
require(grid)
library(psych)
library(knitr)

theme_af <- function(base_size = 15) {
  theme(
    #axis.text = element_text(size = 5),
    #axis.title = element_text(size = 9)
    #axis.line =         theme_blank(),
    axis.text.x =       element_text(size = base_size * 0.8 , lineheight = 0.9, colour = "grey50", angle = 90, vjust=0.5, hjust = 1),
    axis.text.y =       element_text(size = base_size * 0.8, lineheight = 0.9, colour = "grey50", hjust = 1),
    axis.ticks =        element_line(colour = "grey50"),
    axis.title.x =      element_text(size = base_size, vjust = 0.5),
    axis.title.y =      element_text(size = base_size, angle = 90, vjust = 0.5),
    axis.ticks.length = unit(0.15, "cm"),
    
    
    #legend.background = theme_rect(colour="white"), 
    #legend.key =        theme_rect(fill = "grey95", colour = "white"),
    legend.key.size =   unit(1.2, "lines"),
    legend.text =       element_text(size = base_size * 0.8),
    legend.title =      element_text(size = base_size * 0.8, face = "bold", hjust = 0),
    legend.position =   "right",
    
    #panel.background =  theme_rect(fill = "grey90", colour = NA), 
    #panel.border =      theme_blank(), 
    #panel.grid.major =  theme_line(colour = "white"),
    #panel.grid.minor =  theme_line(colour = "grey95", size = 0.25),
    #panel.margin =      unit(0.25, "lines"),
    
    #strip.background =  theme_rect(fill = "grey80", colour = NA), 
    #strip.text.x =      element_text(size = base_size * 0.8),
    #strip.text.y =      element_text(size = base_size * 0.8, angle = -90),
    
    plot.title =        element_text(size = base_size * 1.2)
    #plot.margin =       unit(c(1, 1, 0.5, 0.5), "lines")
    
  )
}
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

# Read stats data
qc_filename="HSC.4vs4.STAR_stats.tsv"

qc=read.table(qc_filename,header = T,stringsAsFactors = T)
rownames(qc)<-qc[,1]
# Summary statistics table
summary(qc)
# Summary by condition
qc_WT = qc[qc$Condition=="WT",]
qc_HET = qc[qc$Condition=="HET",]
# Write file
write.table(describe(qc_WT), "Mapping_statistics/WT.stats", sep="\t", quote=F)
write.table(describe(qc_HET), "Mapping_statistics/HET.stats", sep="\t", quote=F)
write.table(describe(qc), "Mapping_statistics/Total.stats", sep="\t", quote=F)


cols<-brewer.pal(length(c(13:18)), "Pastel1")

### Plot 1 Uniquely mapped reads
title="Total Mapped reads"
yaxis_txt="Mapped reads %"
xaxis_txt="Sample name"

p=ggplot(qc,aes(y=Mapped_reads_perc, x=Sample_name ,fill=Condition))
plot1=p+theme_af()+
  labs(x=xaxis_txt,y=yaxis_txt)+
  ggtitle(title)+
  geom_bar(stat = "identity", alpha = 0.6)

### Plot 2 Uniquely mapped reads
title="Uniquely mapped reads"
yaxis_txt="Uniquely mapped reads %"
xaxis_txt="Sample name"

p=ggplot(qc,aes(y=Uniquely_mapped_reads_perc, x=Sample_name ,fill=Condition))
plot2=p+theme_af()+
  labs(x=xaxis_txt,y=yaxis_txt)+
  ggtitle(title)+
  geom_bar(stat = "identity", alpha = 0.6)

### Plot 3 Multimapping reads
title="Multimapping reads"
yaxis_txt="Multimapping reads %"
xaxis_txt="Sample name"

p=ggplot(qc,aes(y=MultiMapping_reads_perc, x=Sample_name ,fill=Condition))
plot3 = p+theme_af()+
  labs(x=xaxis_txt,y=yaxis_txt)+
  ggtitle(title)+
  geom_bar(stat = "identity", alpha = 0.6)

### Plot 4 Unmapped reads
title="Unmapped reads"
yaxis_txt="Unmapped reads %"
xaxis_txt="Sample name"

p=ggplot(qc,aes(y=Unmapped_reads_perc, x=Sample_name ,fill=Condition))
plot4 = p+theme_af()+
  labs(x=xaxis_txt,y=yaxis_txt)+
  ggtitle(title)+
  geom_bar(stat = "identity", alpha = 0.6)

### Make a multi plot of the statistics
multiplot(plot1,plot2,plot4,plot3,cols=2)
#grid.arrange(plot1, plot2, plot4, plot3, ncol=2)

