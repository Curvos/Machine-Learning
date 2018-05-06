## 载入svm函数包
library('e1071')
## 读取数据
data <- read.table("eset_gcrma.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# 留一法测试函数
LOO <- function(fit.data,ori.data){
  pre <- NULL
  for(i in 1:nrow(ori.data))
  {
    # 逐行选取测试集
    test <- ori.data[i,]
    # 逐行预测
    pre <- c(pre,as.numeric(as.character(predict(fit.data,ori.data[i,]))))
  }
  return(pre)
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
# 0代表对照组Ba，1代表实验组Le
pca_data <- cbind(pca_data,c(rep(0,58),rep(1,58)))
colnames(pca_data) <- c(colnames(pca_data[,-11]),'Class')
pca_data <- as.data.frame(pca_data)

# 使用radial kernel进行SVM分类建模，绘制PC3~PC4的分类图
sv<-svm(Class~.,data=pca_data,cross=5,type='C-classification',kernel='radial')
summary(sv)
png('SVM_plot.png')
plot(sv,pca_data,PC3~PC4)
dev.off()

# 留一法计算性能指标
re = LOO(sv,pca_data)
TN = length(which(re[1:58] == 0)) / 58
FN = length(which(re[1:58] == 1)) / 58
TP = length(which(re[59:116] == 1)) / 58
FP = length(which(re[59:116] == 0)) / 58
TN
FN
TP
FP
