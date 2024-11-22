library(ggplot2)
args = commandArgs(trailingOnly=TRUE)
infile=paste("output",args[1],"relative_affinity_xyplot.input",sep="/")
outfile=paste("output",args[1],"relative_affinity.png",sep="/")
dat <- read.table(infile,header = TRUE)
p<-ggplot(data=dat)+geom_hline(aes(yintercept=0),size=0.2)+geom_vline(aes(xintercept=0),size=0.2)+geom_point(mapping=aes(x=HT_SELEX,y=CAP_SELEX,colour=Class,shape=Class,size=Class))+scale_shape_manual(values=c(21, 16, 16, 16))+labs(title=args[1])+scale_color_manual(values=c('#8B0000','#F8766D','#00BA38','#619CFF'))+scale_size_manual(values=c(1.5, 1, 1, 1))+geom_abline(slope=1,intercept=0, colour="#696969", linetype="dashed")+scale_x_continuous(limit=c(-0.01,1.01),breaks=c(0,0.25,0.5,0.75,1))+scale_y_continuous(limit=c(-0.01,1.01),breaks=c(0,0.25,0.5,0.75,1))+ coord_fixed(ratio=1) + theme(
  panel.background = element_rect(fill = "transparent",colour = NA),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  plot.background = element_rect(fill = "transparent",colour = NA),
#  axis.line.x=element_line(linetype=1,color="black",size=0.2),
#  axis.line.y=element_line(linetype=1,color="black",size=0.2)
  axis.ticks.x=element_line(size=0),
  axis.ticks.y=element_line(size=0),
  plot.title = element_text(hjust = 0.5)
  #axis.text.y = element_text(hjust = 40),
  #axis.text.x = element_text(vjust = 5)
  )
ver <- data.frame(x=c(0,0,0,0,0,0,0,1),y=c(0,0,0,0,0,0,0,1))
p + geom_segment(data=ver,aes(x=c(0,0,0,0,1,0.75,0.5,0.25),y=c(1,0.75,0.5,0.25,0,0,0,0),xend=c(-0.01,-0.01,-0.01,-0.01,1,0.75,0.5,0.25),yend=c(1,0.75,0.5,0.25,-0.01,-0.01,-0.01,-0.01)),size=0.2)
#p + labs(title = args[1])
png(file=outfile,width=2000,height=2000,res=300)
print(p)
dev.off()

