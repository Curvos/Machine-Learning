data <- read.table("last.txt",header=T,sep="\t")

clusting(data[,1:180],"clustering_tree_col.png",3)
clusting(t(data[,1:180]),"clustering_tree_row.png",2)

target <- scale(t(data[,1:180]),center=T,scale=T)
pca <- princomp(target,cor=T)

png("bar-stone_plot1.png",width=600*3,height=3*300,res=72*3) 
par(mfrow=c(1,2),las=2)
#条形图
plot(pca)
abline(h=1,type="2",col="red")
#主成分的碎石图
screeplot(pca, type="lines")
abline(h=1,type="2",col="red")
dev.off()

pca_data <- predict(pca)
clusting(pca_data[,1:4],"clustering_tree_pca.png",8)



data2 <- read.table("fit.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
clusting(data2[,1:180],"PP_cluster_col.png",2)
clusting(t(data2[,1:180]),"PP_cluster_row.png",2)

target1 <- scale(t(data2[,1:180]),center=T,scale=T)
pca1 <- princomp(target1,cor=T)
pca_data1 <- predict(pca1)

png("bar-stone_plot2.png",width=600*3,height=3*300,res=72*3) 
par(mfrow=c(1,2),las=2)
#条形图
plot(pca1)
abline(h=1,type="2",col="red")
#主成分的碎石图
screeplot(pca1, type="lines")
abline(h=1,type="2",col="red")
dev.off()

clusting(pca_data1[,1],"clustering_tree_pca1.png",4)

# 聚类函数
clusting <- function(data,picname,count) {
  d <- dist(data, method = "euclidean")
  hc <- hclust(d,"single")
  png(picname,width=600,height=300)
  plot(hc,xlab="clustering")
  rect.hclust(hc,k=count)
  dev.off()
}