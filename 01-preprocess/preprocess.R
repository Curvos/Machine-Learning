library(affy)
library(tcltk)
library(simpleaffy)
library(affyPLM)
# 载入色彩包
library(RColorBrewer)


## 数据读取
data.raw <- read.affy(covdesc="phenodata.txt", path="data")
file.names=sampleNames(data.raw)

data.raw.gcrma <- gcrma(data.raw)


## 数据质量控制
# 设置调色板
cols <- brewer.pal(8, "Set1")
# 对标准化之前的探针信号强度做箱线图
png("plots\\pre_std.png")
boxplot(data.raw, col=cols)
dev.off()
# 对标准化之后的探针信号强度做箱线图，需要先安装 affyPLM 包，以便解析 data.raw.gcrma 数据
png("plots\\post_std.png")

boxplot(data.raw.gcrma, col=cols)
dev.off()
# 标准化前后的箱线图会有些变化
# 但是密度曲线图看起来更直观一些
# 对标准化之前的数据做密度曲线图
png("plots\\pre_std_curve.png")
hist(data.raw, col=cols)
dev.off()
# 对标准化之后的数据做密度曲线图
png("plots\\post_std_curve.png")
hist(data.raw.gcrma, col=cols)
dev.off()

# 从 CEL 文件读取探针信号强度:
data.raw.qc <- fitPLM(data.raw)
# RLE (Relative Log Expression 相对表达量取对数) 图中
# 所有的值都应该接近于零。 GSM524665.CEL 芯片数据由于有人为误差而例外
png("plots\\RLE.png")
RLE(data.raw.qc, main="RLE")
dev.off()
# 也可以用 NUSE (Normalised Unscaled Standard Errors)作图比较.
# 对于绝大部分基因，标准差的中位数应该是1。
# 芯片 GSM524665.CEL 在这个图中，同样是一个例外
png("plots\\NUSE.png")
NUSE(data.raw.qc, main="NUSE")
dev.off()

## 灰度图扫描
# 芯片数量
n.cel <- length(file.names)
#par(mfrow = c(ceiling(n.cel/2), 2))
#par(mar = c(0.5, 0.5, 2, 0.5))
# 设置调色板颜色为灰度
pallette.gray <- c(rep(gray(0:10/10), times = seq(1, 41, by = 4)))
# 通过for循环逐个作图
for (i in 1:n.cel){
  out_name = paste("gray_plots\\", file.names[i], sep = "")
  out_name = paste(out_name, ".png", sep = "")
  png(out_name)
  image(data.raw[, i], col = pallette.gray)
  dev.off()
} 


## 背景处理
data.rma <- bg.correct(data.raw, method="rma")
data.mas <- bg.correct(data.raw, method="mas")
class(data.rma)
class(data.mas)
class(data.raw)
head(pm(data.raw)-pm(data.mas),2)
head(pm(data.raw)-pm(data.rma),2)
head(mm(data.raw)-mm(data.mas),2)
head(mm(data.raw)-mm(data.rma),2)
identical(mm(data.raw), mm(data.rma))


## 归一化处理
# loess
data.rma.lo <- normalize(data.rma, method = "loess")
# 非线性缩放方法
# data.mas.nl <- normalize(data.mas, method = "invariantset")

## 汇总
eset.rma.liwong <- computeExprSet(data.rma.lo, pmcorrect.method="pmonly",
                                  summary.method="liwong")


out = eset.rma.liwong@assayData$exprs
write.table(out, file = "result.txt", quote = F, sep = '\t')