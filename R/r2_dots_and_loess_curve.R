#######
####### Output a graph with GGplot2 with r2 as a function of the distance between markers (LD calculation) 
#######

#plot r2 
setwd(".") #the file has to be in the current directory
compa<-read.table("LD_Chrall_except12_Haploview.reduced.txt",sep="\t",header=T)
colnames(compa)<-c("l1","l2","dprime", "lod", "r2", "cilow","cihi","dist","tint")

library(ggplot2)
plotr2<-ggplot(data=compa, aes(x=dist, y=r2)) + geom_point() + stat_smooth(method="auto", fill = "grey50", size = 2, alpha = 1) + # pick a binwidth that is not too small 
#facet_wrap(~ chr,ncol=2,scales="free_x") + # seperate plots for each chr, x-scales can differ from chr to chr
ggtitle("intra chromosomal LD (r2) decay of marker pairs over all chromosomes as a function of the distance between pairs") +
#ylim(0,5000) +
#xlim(0,5) +
xlab("distance") + 
ylab("r2") + 
theme_bw() # black and white theme

# save the plot to .png file
png("r2_plots_and_loess.png",3000,1500)
print(plotr2)
dev.off()

#plot(plotr2)



