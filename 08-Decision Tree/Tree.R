## 载入svm和knn函数包
library('e1071')
library('caret')
library('class')
library('partykit')
library('randomForest')
library('rpart')
library('discretization')

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


## svm-rfe 特征选择
df.data <- as.data.frame(data.set)
d.train <- as.data.frame(train)
set.seed(123)
ctrl= rfeControl(functions = rfFuncs, method = "cv",verbose = FALSE, 
                 returnResamp = "final")
svm.features2 <- rfe(data.set[,1:421], data.set[, 422],
                    sizes = c(50,100,150,200,250,300,350,400),
                    rfeControl = ctrl)
#rfeCNTL <- rfeControl(functions = lrFuncs, method = "cv", number= 10)
#svm.features <- rfe(d.train[,1:421], d.train[, 422],
#                    sizes = c(50,100,150,200,250,300,350,400),
#                    rfeControl = ctrl)
#method = "svmLinear")

print(svm.features2)
sel.cols <- svm.features2$optVariables
ranked.svm <- subset(df.data, select = sel.cols)
ranked.svm$Tag <- df.data$Tag

acc.svm.t <- NULL
acc.kmm.t <- NULL
for (i in 20:116){
  n.samples = sample(1:116,i)
  n.test = ranked.svm[n.samples,]
  sv.f <- svm(Tag~., data=n.test, cross=5, type='C-classification', kernel='radial')
  acc.svm.t = c(acc.svm.t, sv.f$tot.accuracy)
  
  ind.f <- sample(2, nrow(n.test), replace=T, prob=c(0.7,0.3))
  train.f <- n.test[ind.f == 1, ]
  test.f <- n.test[ind.f == 2, ]
  knn.test.f <- knn(train.f[,-101], test.f[,-101], train.f[,101], k=9)
  tab <- table(knn.test.f, test.f[,101])
  acc.kmm.t <- c(acc.kmm.t,(tab[1,1]+tab[2,2]) / nrow(test.f))
}
# SVM plot
png('SVM after rfe rank.png')
plot(x=seq(20,116), y=acc.svm.t, xlab='n from 20 to 116', ylab="Accuracy", main='SVM after rfe rank', 
     col='olivedrab')
dev.off()
# KNN plot
png('KNN after rfe rank.png')
plot(x=seq(20,116), y=acc.kmm.t, xlab='n from 20 to 116', ylab="Accuracy", main='KNN after rfe rank', 
     col='orangered')
dev.off()

# 对比svm-Rfe及主成分分析模型
wilcox.test(acc.kmm.pca, acc.kmm.t, paired = TRUE, alternative = "less")
wilcox.test(acc.svm.pca, acc.svm.t, paired = TRUE, alternative = "less")


## 决策树
# 决策树测试
set.seed(123)
dis_result <- mdlp(n.test)
dis.df = dis_result$Disc.data
tree.biop2 <- rpart(Tag ~ ., data = dis.df)
tpp <- tree.biop2$cptable
row.names(tpp) <- NULL
best.rel <- which.min(tpp[,3])
cp <- min(tree.biop2$cptable[best.rel, ])
prune.tree.biop2 <- prune(tree.biop2, cp = cp)
plot(as.party(prune.tree.biop2))

# 对整张表达谱进行决策树建立
dis.result <- mdlp(df.data)
dis.df <- dis_result$Disc.data
dis.df$Tag <- factor(dis.df$Tag)

acc.tree <- NULL
for (i in 20:116){
  n.samples = sample(1:116,i)
  n.test = dis.df[n.samples,]
  ind.f <- sample(2, nrow(n.test), replace=T, prob=c(0.7,0.3))
  biop.train <- n.test[ind.f == 1, ]
  biop.test <- n.test[ind.f == 2, ]
  
  tree.biop <- rpart(Tag ~ ., data = biop.train)
  tpp <- tree.biop$cptable
  row.names(tpp) <- NULL
  best.rel <- which.min(tpp[,3])
  cp <- min(tpp[best.rel,])
  prune.tree.biop <- prune(tree.biop, cp = cp)
  rparty.test <- predict(prune.tree.biop, newdata = biop.test, type
                         = "class")
  tab <- table(rparty.test, biop.test$Tag)
  # 计算准确率
  acc.tree <- c(acc.tree,(tab[1,1]+tab[2,2]) / nrow(biop.test))
}
# Decision Tree plot
png('Decision Tree.png')
plot(x=seq(20,116), y=acc.tree, xlab='n from 20 to 116', ylab="Accuracy", 
     main='Decision tree after discretization', col='dodgerblue')
dev.off()

# 以决策树模型对比KNN和SVM模型(pca后)
wilcox.test(acc.kmm.pca, acc.tree, paired = TRUE)
wilcox.test(acc.svm.pca, acc.tree, paired = TRUE)
