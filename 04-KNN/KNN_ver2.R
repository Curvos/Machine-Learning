## 读取数据
data <- read.table("eset_gcrma.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

## 编写所需函数
# clust：分别接收数据集、训练集和K值
clust <- function(dataset,testset,K) {
  distance <- NULL
  for (i in 1:nrow(dataset)) {
    # 计算距离
    distance <- c(distance,dist(rbind(testset,dataset[i,]),method = "euclidean"))
  }
  names(distance) <- rownames(dataset)
  temp <- sort(distance)
  result <- temp[1:10]
  ba <- 0
  le <- 0
  for(i in 1:K) {
    if(names(result[i])=="Ba") {
      ba = ba + 1
    }
    else {
      le = le + 1
    }
  }
  if (ba >= le) {
    return("Ba")
  }
  else  {
    return("Le")
  }
}

# 循环实现KNN算法
KNN <- function(count,K) {
  KNN_resu <- NULL
  for(i in 1:count) {
    temp <- NULL
    # 随机取一行样本作为测试集
    t <- sample(1:116, 1)
    if(t >= 1 & t <= 58) {
      temp <- "Ba"
    }
    else{
      temp <- "Le"
    }
    res <- clust(pca_data[-t,],pca_data[t,],K)
    KNN_resu <- c(KNN_resu,paste(temp,res,sep="-"))
  }
  # 敏感度
  TR <- length(which(KNN_resu=="Ba-Ba" | KNN_resu=="Le-Le" ))
  accuracy <- TR/count
  # Ba检测准确率
  Ba_true <- length(which(KNN_resu=="Ba-Ba"))
  Ba_all <- length(which(KNN_resu=="Ba-Ba" | KNN_resu=="Le-Ba" ))
  Ba_TRP <- Ba_true/Ba_all
  # Le检测准确率
  Le_ture <- length(which(KNN_resu=="Le-Le"))
  Le_all <- length(which(KNN_resu=="Le-Le" | KNN_resu=="Ba-Le"))
  Le_TRP <- Le_ture/Le_all
  
  result <- c(accuracy,Ba_TRP,Le_TRP)
  names(result) <- c("accuracy","Ba_TRP","Le_TRP")
  
  return(result)
}


## 数据降维，为 Ba：正常组和 Le：来曲唑组贴上标签
# 进行主成分分析，绘制滚石图
pca <- prcomp(t(data),scale=FALSE)
png("bar_stone.png") 
par(mfrow=c(1,2),las=2)
plot(pca)
abline(h=1,type="2",col="red")
screeplot(pca, type="lines")
abline(h=1,type="2",col="red")
dev.off()
pca_data <- predict(pca)[,1:10]
#  Ba：正常组 Le：来曲唑予药组
rownames(pca_data) <- c(rep("Ba",58),rep("Le",58))



## 使用KNN算法对116个样本进行分类，K取9
aim <- KNN(116,9)
for(i in 1:10) {
  aim <- rbind(aim,KNN(116,9))
}
row.names(aim) = NULL
result = as.data.frame(apply(aim[,1:3],2,mean))
write.table(result, file = 'result.txt')

