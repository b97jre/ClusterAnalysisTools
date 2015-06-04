
library(ggplot2)

chr1 = 30427671
chr2 = 19698289
chr3 = 23459830
chr4 = 18585056
chr5 = 26975502

fileName = '/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/Wellmer2007Clusters.txt'
ClusterInfo = read.table(file = fileName,header = TRUE)
ClusterInfo$Source = "SEG"

fileName = '/Users/johanreimegard/Vetenskap/Data/clusterProject/TAIR10/Wellmer2004Clusters.txt'
ClusterInfo2 = read.table(file = fileName,header = TRUE)
ClusterInfo2$Source = "MS1"

ClusterInfo = rbind(ClusterInfo, ClusterInfo2)

ClusterInfo$Location = 1
ClusterInfo[ClusterInfo$Chr == 'Chr1', 'Location'] = ClusterInfo[ClusterInfo$Chr == 'Chr1', 'Start']/chr1 
ClusterInfo[ClusterInfo$Chr == 'Chr2', 'Location'] = ClusterInfo[ClusterInfo$Chr == 'Chr2', 'Start']/chr2 
ClusterInfo[ClusterInfo$Chr == 'Chr3', 'Location'] = ClusterInfo[ClusterInfo$Chr == 'Chr3', 'Start']/chr3 
ClusterInfo[ClusterInfo$Chr == 'Chr4', 'Location'] = ClusterInfo[ClusterInfo$Chr == 'Chr4', 'Start']/chr4 
ClusterInfo[ClusterInfo$Chr == 'Chr5', 'Location'] = ClusterInfo[ClusterInfo$Chr == 'Chr5', 'Start']/chr5 



ChromosomePlot = ggplot(data = NULL)
ChromosomePlot + 
  scale_y_discrete(name="",  limits = rev(c("Chr1","Chr2","Chr3","Chr4","Chr5"))) +
  scale_x_continuous(limits=c(0,1), breaks=NULL, name="") +
  geom_point(data=ClusterInfo, mapping=aes(y=Chr, x=Location , color = Source,size=ClusterInfo$total), alpha  = 0.66)+
  labs(title = "Cluster positions")

ggsave("ClusterDistributions.pdf",width = 4, height = 3)




