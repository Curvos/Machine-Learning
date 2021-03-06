## 载入svm和knn函数包
library('e1071')
library('caret')
library('class')


# 留一法测试函数
LOO <- function(fit.data,ori.data){
  pre <- NULL
  for(i in 1:nrow(ori.data))
  {
    # 逐行选取测试集
    test_LOO <- ori.data[i,]
    # 逐行预测
    pre <- c(pre,as.numeric(as.character(predict(fit.data,ori.data[i,]))))
  }
  return(pre)
}

## 差异表达基因数据读取
data <- read.table("t_test.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
raw.profile <- read.table("eset_rma.txt", header = TRUE, sep = "\t", 
                          stringsAsFactors = FALSE)
# 对差异表达谱以fold.change绝对值从高到低进行排序
expr_sorted <- data[order(abs(data[,3]),decreasing=T),]
# 在整张表达谱上选取差异表达基因
dataset <- raw.profile[row.names(expr_sorted),]
# 对数据进行转置
data.set <- t(dataset)
# 检测选取后是否依旧符合排序
row.names(dataset) == row.names(expr_sorted)
# 添加tag, 0代表对照组Ba，1代表实验组Le
data.set <- cbind(data.set,c(rep(0,58),rep(1,58)))
colnames(data.set) <- c(colnames(data.set[,-422]),"Tag")


## BEFORE pca
# svm 测试
sv <- svm(Tag~., data=data.set, cross=5, type='C-classification', kernel='radial')
summary(sv)
# knn

# 训练集：测试集 = 3：1 测试
ind <- sample(2, nrow(data.set), replace=T, prob=c(0.7,0.3))
train <- data.set[ind == 1, ]
test <- data.set[ind == 2, ]
knn.test <- knn(train[,-422], test[,-422], train[,422], k=9)


acc.svm <- NULL
acc.kmm <- NULL
for (i in 20:116){
  #n.samples = sample(data.set, 100, replace = F)
  #n.samples = sample(1:116,20)
  #n.test = data.set[n.samples,]
  n.samples = sample(1:116,i)
  n.test = data.set[n.samples,]
  sv.f <- svm(Tag~., data=n.test, cross=5, type='C-classification', kernel='radial')
  acc.svm = c(acc.svm, sv.f$tot.accuracy)
  
  ind.f <- sample(2, nrow(n.test), replace=T, prob=c(0.7,0.3))
  train.f <- n.test[ind.f == 1, ]
  test.f <- n.test[ind.f == 2, ]
  knn.test.f <- knn(train.f[,-422], test.f[,-422], train.f[,422], k=9)
  tab <- table(knn.test.f, test.f[,422])
  acc.kmm <- c(acc.kmm,(tab[1,1]+tab[2,2]) / nrow(test.f))
}
# SVM plot
png('SVM before pca.png')
plot(x=seq(20,116), y=acc.svm, xlab='n from 20 to 116', ylab="Accuracy", main='SVM before pca', 
     col='olivedrab')
dev.off()
# KNN plot
png('KNN before pca.png')
plot(x=seq(20,116), y=acc.kmm, xlab='n from 20 to 116', ylab="Accuracy", main='KNN before pca', 
     col='orangered')
dev.off()

## AFTER pca
# 进行主成分分析，绘制滚石图
pca <- prcomp(data.set,scale=FALSE)
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
colnames(pca_data) <- c(colnames(pca_data[,-11]),'Tag')
pca_data <- as.data.frame(pca_data)

acc.svm.pca <- NULL
acc.kmm.pca <- NULL
for (i in 20:116){
  #n.samples = sample(data.set, 100, replace = F)
  #n.samples = sample(1:116,20)
  #n.test = data.set[n.samples,]
  n.samples = sample(1:116,i)
  n.test = pca_data[n.samples,]
  sv.f <- svm(Tag~., data=n.test, cross=5, type='C-classification', kernel='radial')
  acc.svm.pca = c(acc.svm.pca, sv.f$tot.accuracy)
  
  ind.f <- sample(2, nrow(n.test), replace=T, prob=c(0.7,0.3))
  train.f <- n.test[ind.f == 1, ]
  test.f <- n.test[ind.f == 2, ]
  knn.test.f <- knn(train.f[,-11], test.f[,-11], train.f[,11], k=9)
  tab <- table(knn.test.f, test.f[,11])
  acc.kmm.pca <- c(acc.kmm.pca,(tab[1,1]+tab[2,2]) / nrow(test.f))
}

# SVM.pca plot
png('SVM after pca.png')
plot(x=seq(20,116), y=acc.svm.pca, xlab='n from 20 to 116', ylab="Accuracy", main='SVM after pca', 
     col='olivedrab')
dev.off()
# KNN.pca plot
png('KNN after pca.png')
plot(x=seq(20,116), y=acc.kmm.pca, xlab='n from 20 to 116', ylab="Accuracy", main='KNN after pca', 
     col='orangered')
dev.off()