#install.packages("glmnet")
#install.packages("pROC")

#install.packages("glmnet")
#install.packages("pROC")
#install.packages("ggplot2")

library(glmnet)
library(pROC)
library(ggplot2)

# =======================
# 输入区：你只改这里
# =======================
expFile      = "train.normalize.txt"          # 内部/训练表达矩阵：gene x sample
expFile_ext  = "external.normalize.txt"       # 外部验证表达矩阵：gene x sample；不需要外部就留空 "" 或 NULL
geneFile     = "interGenes.txt"               # 模型基因列表（每行一个基因）
setwd("C:\\Users\\1453621405\\Desktop\\MDD\\final\\ROC")

set.seed(123)

# =======================
# 工具函数：在“定义阈值”下计算临床指标（TP/TN/FP/FN, Se/Sp/PPV/NPV等）
# =======================
calc_clinical_metrics <- function(y_true, prob, thr){
  pred_class <- ifelse(prob >= thr, 1, 0)
  
  TP <- sum(pred_class==1 & y_true==1)
  TN <- sum(pred_class==0 & y_true==0)
  FP <- sum(pred_class==1 & y_true==0)
  FN <- sum(pred_class==0 & y_true==1)
  
  sens <- ifelse((TP+FN)==0, NA, TP/(TP+FN))
  spec <- ifelse((TN+FP)==0, NA, TN/(TN+FP))
  ppv  <- ifelse((TP+FP)==0, NA, TP/(TP+FP))
  npv  <- ifelse((TN+FN)==0, NA, TN/(TN+FN))
  acc  <- ifelse((TP+TN+FP+FN)==0, NA, (TP+TN)/(TP+TN+FP+FN))
  
  data.frame(
    TP=TP, TN=TN, FP=FP, FN=FN,
    Sensitivity=sens,
    Specificity=spec,
    PPV=ppv,
    NPV=npv,
    Accuracy=acc
  )
}

# =======================
# 读取内部数据
# =======================
rt = read.table(expFile, header=TRUE, sep="\t", check.names=FALSE, row.names=1)

# 分组标签
y = gsub("(.*)\\_(.*)", "\\2", colnames(rt))
y = ifelse(y=="Control", 0, 1)

# 读基因列表
geneRT = read.table(geneFile, header=FALSE, sep="\t", check.names=FALSE)

# =======================
# 单基因ROC（原样保留）
# =======================
bioCol = rainbow(nrow(geneRT), s=0.9, v=0.9)
aucText = c()
k = 0
for(x in as.vector(geneRT[,1])){
  k = k + 1
  roc1 = roc(y, as.numeric(rt[x, ]))
  if(k==1){
    pdf(file="ROC.genes.pdf", width=5, height=4.75)
    plot(roc1, print.auc=FALSE, col=bioCol[k], legacy.axes=TRUE, main="")
    aucText = c(aucText, paste0(x,", AUC=",sprintf("%.3f", roc1$auc[1])))
  }else{
    plot(roc1, print.auc=FALSE, col=bioCol[k], legacy.axes=TRUE, main="", add=TRUE)
    aucText = c(aucText, paste0(x,", AUC=",sprintf("%.3f", roc1$auc[1])))
  }
}
legend("bottomright", aucText, lwd=2, bty="n", col=bioCol[1:nrow(geneRT)])
dev.off()

# =======================
# 逻辑回归模型（原主体保留 + glm predict修正）
# =======================
rt_m = rt[as.vector(geneRT[,1]), , drop=FALSE]
rt_m = as.data.frame(t(rt_m))
rt_m$y = y

logit = glm(y ~ ., family=binomial(link='logit'), data=rt_m)

# 内部预测概率
pred = predict(logit, newdata=rt_m, type="response")

# 内部ROC（AUC+CI）
roc_model = roc(y, as.numeric(pred))
ci1 = ci.auc(roc_model, method="bootstrap", boot.n=10000)
ciVec = as.numeric(ci1)

# =======================
# 定义阈值：训练集Youden阈值（defined threshold）
# =======================
thr = as.numeric(coords(roc_model, x="best", best.method="youden",
                        ret="threshold", transpose=FALSE))

# 在定义阈值下计算临床指标（内部）
clin_int = calc_clinical_metrics(y_true=y, prob=as.numeric(pred), thr=thr)

# 画内部模型ROC
pdf(file="ROC.model.pdf", width=5, height=4.75)
plot(roc_model, print.auc=TRUE, col="red", legacy.axes=TRUE, main="Model (Internal)")
text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec[1]),"-",sprintf("%.03f",ciVec[3])), col="red")
text(0.39, 0.36, paste0("Defined threshold (Youden): ", sprintf("%.6f", thr)), col="red")
dev.off()

# =======================
# 输出：审稿人要求的临床可解释指标（Internal/External）
# =======================
metrics_clinical = data.frame(
  Cohort = "Internal",
  AUC = as.numeric(auc(roc_model)),
  CI_low = ciVec[1],
  CI_high = ciVec[3],
  Threshold = thr,
  clin_int
)

# =======================
# 外部验证：阈值不重选，固定用训练阈值thr
# =======================
if(!is.null(expFile_ext) && expFile_ext != "" && file.exists(expFile_ext)){
  rt2 = read.table(expFile_ext, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
  y2 = gsub("(.*)\\_(.*)", "\\2", colnames(rt2))
  y2 = ifelse(y2=="Control", 0, 1)
  
  rt2_m = rt2[as.vector(geneRT[,1]), , drop=FALSE]
  rt2_m = as.data.frame(t(rt2_m))
  rt2_m$y = y2
  
  pred2 = predict(logit, newdata=rt2_m, type="response")
  
  roc_ext = roc(y2, as.numeric(pred2))
  ci2 = ci.auc(roc_ext, method="bootstrap", boot.n=10000)
  ciVec2 = as.numeric(ci2)
  
  clin_ext = calc_clinical_metrics(y_true=y2, prob=as.numeric(pred2), thr=thr)
  
  metrics_clinical = rbind(metrics_clinical, data.frame(
    Cohort = "External",
    AUC = as.numeric(auc(roc_ext)),
    CI_low = ciVec2[1],
    CI_high = ciVec2[3],
    Threshold = thr,
    clin_ext
  ))
  
  # 外部ROC图
  pdf(file="ROC.model.external.pdf", width=5, height=4.75)
  plot(roc_ext, print.auc=TRUE, col="blue", legacy.axes=TRUE, main="Model (External)")
  text(0.39, 0.43, paste0("95% CI: ",sprintf("%.03f",ciVec2[1]),"-",sprintf("%.03f",ciVec2[3])), col="blue")
  text(0.39, 0.36, paste0("Defined threshold (train Youden): ", sprintf("%.6f", thr)), col="blue")
  dev.off()
  
  # 内外叠加ROC
  pdf(file="ROC.model.internal_vs_external.pdf", width=5.5, height=4.75)
  plot(roc_model, col="red", legacy.axes=TRUE, lwd=2, main="Model ROC: Internal vs External")
  plot(roc_ext, col="blue", lwd=2, add=TRUE)
  legend("bottomright",
         legend=c(paste0("Internal AUC=", sprintf("%.3f", as.numeric(auc(roc_model)))),
                  paste0("External AUC=", sprintf("%.3f", as.numeric(auc(roc_ext))))),
         col=c("red","blue"), lwd=2, bty="n")
  dev.off()
}

# 输出csv（临床指标完整版）
write.csv(metrics_clinical, file="ROC.metrics.clinical.csv", row.names=FALSE)
print(metrics_clinical)

# =======================
# 新增：柱状图（Internal vs External 的 Se/Sp/PPV/NPV）
# =======================
plot_metrics <- c("Sensitivity","Specificity","PPV","NPV")

bar_df <- metrics_clinical[, c("Cohort", plot_metrics)]
bar_long <- data.frame(
  Cohort = rep(bar_df$Cohort, each=length(plot_metrics)),
  Metric = rep(plot_metrics, times=nrow(bar_df)),
  Value  = as.vector(t(as.matrix(bar_df[, plot_metrics])))
)

# 去掉 NA（比如 NPV 分母为0时会NA）
bar_long <- bar_long[!is.na(bar_long$Value), ]

pdf("Clinical_metrics.barplot.pdf", width=7, height=4.5)
ggplot(bar_long, aes(x=Metric, y=Value, fill=Cohort)) +
  geom_col(position="dodge", width=0.7) +
  coord_cartesian(ylim=c(0,1)) +
  labs(title=paste0("Clinical metrics at defined threshold (", sprintf("%.4f", thr), ")"),
       x="", y="Value") +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))
dev.off()

