library(affy)
library(tcltk)
library(simpleaffy)
library(affyPLM)

## 数据读取
data.raw <- read.affy(covdesc="phenodata.txt", path="data")
file.names=sampleNames(data.raw)
n.cel <- length(data.raw)
data.raw.gcrma <- gcrma(data.raw)


## 使用rma方法进行预处理
eset.rma <- rma(data.raw)


## 计算基因表达量
# 计算平均值，并做对数转换
emat.rma.log2 <- exprs(eset.rma)
emat.rma.nologs <- 2^emat.rma.log2
class(emat.rma.log2)
head(emat.rma.log2, 1)

# 使用rma分别计算各组平均表达量和差异表达倍数
# 第1：20列为 ISb 组，与21：40列 ISs 组互为对照
# 第41：60列为 IRb 组，与61：80列 IRs 组互为对照
# 第81：95列为 DBb 组，与96：110列 DBs 组互为对照


# 计算表达均值并经行log转换
# ISb组
results.rma.ISb<- data.frame(apply(emat.rma.log2[,c(1:20)],1,mean))
colnames(results.rma.ISb) = c("ISb")
# ISs组
results.rma.ISs<- data.frame(apply(emat.rma.log2[,c(21:40)],1,mean))
colnames(results.rma.ISs) = c("ISs")
# IRb组
results.rma.IRb<- data.frame(apply(emat.rma.log2[,c(41:60)],1,mean))
colnames(results.rma.IRb) = c("IRb")
# IRs组
results.rma.IRs<- data.frame(apply(emat.rma.log2[,c(61:80)],1,mean))
colnames(results.rma.IRs) = c("IRs")
# DBb组
results.rma.DBb<- data.frame(apply(emat.rma.log2[,c(81:95)],1,mean))
colnames(results.rma.DBb) = c("DBb")
# DBs组
results.rma.DBs<- data.frame(apply(emat.rma.log2[,c(96:110)],1,mean))
colnames(results.rma.DBs) = c("DBs")
# 数据集合并
results.rma = cbind(results.rma.ISb,results.rma.ISs,results.rma.IRb,results.rma.IRs,
                    results.rma.DBb,results.rma.DBs)

# 计算差异倍数
results.rma$fc.IS <- results.rma[,2]-results.rma[,1]
results.rma$fc.IR <- results.rma[,4]-results.rma[,2]
results.rma$fc.DB <- results.rma[,6]-results.rma[,5]
head(results.rma, 2)


## 选取表达基因
data.mas5calls <- mas5calls(data.raw)
# 继续用exprs计算“表达”量
eset.mas5calls <- exprs(data.mas5calls)
head(eset.mas5calls)

# 基因筛选
AP <- apply(eset.mas5calls, 1, function(x)any(x=="P"))
present.probes <- names(AP[AP])
paste(length(present.probes),"/",length(AP))
# 删除中间数据
rm(data.mas5calls)
rm(eset.mas5calls)
# 数据子集提取
results.present <- results.rma[present.probes,]


## 获取差异表达基因

# 计算表达变化超过2倍的基因
sum(abs(results.present[,"fc.IS"])>=0.5)
sum(abs(results.present[,"fc.IR"])>=0.5)
sum(abs(results.present[,"fc.DB"])>=0.5)

## IS 组
# T 检验
apply(abs(results.present[,7:9]), 2, max)
# 选取这些基因
results.st.IS <- results.present[abs(results.present$fc.IS)>=1,]
sel.genes.IS <- row.names(results.st.IS)
# t测验，并选出需要的差异表达基因：
p.value <- apply(emat.rma.log2[sel.genes.IS,], 1, function(x){t.test(x[1:20], x[21:40])$p.value})
results.st.IS$p.value <- p.value
names(results.st.IS)
results.st.IS <- results.st.IS[, c(1,2,7,10)]
# 分别筛选出p<0.05,p<0.01,p<0.001的差异表达基因
results.st.IS.p05 <- results.st.IS[p.value<0.05,]
results.st.IS.p01 <- results.st.IS[p.value<0.01,]
results.st.IS.p001 <- results.st.IS[p.value<0.001,]
head(results.st.IS, 2)
# 查看筛选所得数据数目
nrow(results.st.IS.p05)
nrow(results.st.IS.p01)
nrow(results.st.IS.p001)


## IR 组
# T 检验
apply(abs(results.present[,7:9]), 2, max)
# 选取这些基因
results.st.IR <- results.present[abs(results.present$fc.IR)>=1,]
sel.genes.IR <- row.names(results.st.IR)
# t测验，并选出需要的差异表达基因：
p.value <- apply(emat.rma.log2[sel.genes.IR,], 1, function(x){t.test(x[41:60], x[61:80])$p.value})
results.st.IR$p.value <- p.value
names(results.st.IR)
results.st.IR <- results.st.IR[, c(1,2,7,10)]
# 分别筛选出p<0.05,p<0.01,p<0.001的差异表达基因
results.st.IR.p05 <- results.st.IR[p.value<0.05,]
results.st.IR.p01 <- results.st.IR[p.value<0.01,]
results.st.IR.p001 <- results.st.IR[p.value<0.001,]
head(results.st.IR, 2)
# 查看筛选所得数据数目
nrow(results.st.IR.p05)
nrow(results.st.IR.p01)
nrow(results.st.IR.p001)


## DB 组
# T 检验
apply(abs(results.present[,7:9]), 2, max)
# 计算表达变化超过2倍的基因
# 选取这些基因
results.st.DB <- results.present[abs(results.present$fc.DB)>=1,]
sel.genes.DB <- row.names(results.st.DB)
# t测验，并选出需要的差异表达基因：
p.value <- apply(emat.rma.log2[sel.genes.DB,], 1, function(x){t.test(x[81:95], x[96:110])$p.value})
results.st.DB$p.value <- p.value
names(results.st.DB)
results.st.DB <- results.st.DB[, c(1,2,7,10)]
# 分别筛选出p<0.05,p<0.01,p<0.001的差异表达基因
results.st.DB.p05 <- results.st.DB[p.value<0.05,]
results.st.DB.p01 <- results.st.DB[p.value<0.01,]
results.st.DB.p001 <- results.st.DB[p.value<0.001,]
head(results.st.DB, 2)
# 查看筛选所得数据数目
nrow(results.st.DB.p05)
nrow(results.st.DB.p01)
nrow(results.st.DB.p001)


## 数据输出
write.table(results.st.IS.p05, file = 'output\\IS_p05.txt', sep = '\t')
write.table(results.st.IS.p01, file = 'output\\IS_p01.txt', sep = '\t')
write.table(results.st.IS.p001, file = 'output\\IS_p001.txt', sep = '\t')

write.table(results.st.IR.p05, file = 'output\\IR_p05.txt', sep = '\t')
write.table(results.st.IR.p01, file = 'output\\IR_p01.txt', sep = '\t')
write.table(results.st.IR.p001, file = 'output\\IR_p001.txt', sep = '\t')

write.table(results.st.DB.p05, file = 'output\\DB_p05.txt', sep = '\t')
write.table(results.st.DB.p01, file = 'output\\DB_p01.txt', sep = '\t')
write.table(results.st.DB.p001, file = 'output\\DB_p001.txt', sep = '\t')