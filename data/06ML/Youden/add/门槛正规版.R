library(pROC)
library(ggplot2)

# =======================
# 输入区：只改这里
# =======================
setwd("E:\\AAA_SCI\\Nitrofurantoin_IPF\\06ML\\Youden\\add")

trainFile <- "data.train.txt"
testFile  <- "data.test.txt"

# bootstrap 次数
boot_n_auc <- 10000

# 是否用训练集均值填补缺失
do_impute_by_train_mean <- TRUE

# 是否按训练集做 Z-score 标准化
do_zscore_by_train <- FALSE

# 从样本名提取数据集名称
extract_dataset_name <- function(x){
  sub("_.*$", "", x)
}

set.seed(123)

# =======================
# 读取 train/test（行=样本，列=特征+Type）
# =======================
train_df <- read.table(trainFile, header = TRUE, sep = "\t",
                       check.names = FALSE, row.names = 1)
test_df  <- read.table(testFile,  header = TRUE, sep = "\t",
                       check.names = FALSE, row.names = 1)

if(!("Type" %in% colnames(train_df))) stop("trainFile 缺少 Type 列")
if(!("Type" %in% colnames(test_df)))  stop("testFile 缺少 Type 列")

# =======================
# 分离 Type 和特征
# =======================
y_train <- as.numeric(as.character(train_df$Type))
y_test_all <- as.numeric(as.character(test_df$Type))

if(any(is.na(y_train))) stop("trainFile 的 Type 列无法正确转换为数值，请确认 Type 为 0/1")
if(any(is.na(y_test_all))) stop("testFile 的 Type 列无法正确转换为数值，请确认 Type 为 0/1")

x_train_raw <- train_df[, setdiff(colnames(train_df), "Type"), drop = FALSE]
x_test_raw  <- test_df[,  setdiff(colnames(test_df),  "Type"), drop = FALSE]

x_train_raw[] <- lapply(x_train_raw, function(z) as.numeric(as.character(z)))
x_test_raw[]  <- lapply(x_test_raw,  function(z) as.numeric(as.character(z)))

# 对齐共同特征
common_vars <- intersect(colnames(x_train_raw), colnames(x_test_raw))
if(length(common_vars) == 0) stop("train/test 没有共同特征列")

x_train_raw <- x_train_raw[, common_vars, drop = FALSE]
x_test_raw  <- x_test_raw[,  common_vars, drop = FALSE]

# 去掉全 NA 列
all_na_cols <- colnames(x_train_raw)[
  colSums(!is.na(x_train_raw)) == 0 | colSums(!is.na(x_test_raw)) == 0
]
if(length(all_na_cols) > 0){
  x_train_raw <- x_train_raw[, setdiff(colnames(x_train_raw), all_na_cols), drop = FALSE]
  x_test_raw  <- x_test_raw[,  setdiff(colnames(x_test_raw),  all_na_cols), drop = FALSE]
}

if(ncol(x_train_raw) == 0) stop("过滤后没有可用特征列")

# =======================
# 按数据集名称拆分 test
# =======================
test_sample_names <- rownames(test_df)
test_dataset_names <- extract_dataset_name(test_sample_names)

if(any(is.na(test_dataset_names) | test_dataset_names == "")){
  stop("有样本无法提取数据集名称，请检查 extract_dataset_name() 规则")
}

test_groups <- split(seq_along(test_dataset_names), test_dataset_names)

split_info <- data.frame(
  Sample = test_sample_names,
  Type = y_test_all,
  Dataset = test_dataset_names,
  stringsAsFactors = FALSE
)
write.csv(split_info, "Test_split_by_dataset.csv", row.names = FALSE)

# =======================
# 函数区
# =======================
impute_by_train_mean <- function(x_train, x_new){
  x_train2 <- x_train
  x_new2   <- x_new
  
  train_means <- numeric(ncol(x_train2))
  names(train_means) <- colnames(x_train2)
  
  for(j in seq_len(ncol(x_train2))){
    mu <- mean(as.numeric(x_train2[, j]), na.rm = TRUE)
    if(is.na(mu)) mu <- 0
    train_means[j] <- mu
    
    idx1 <- which(is.na(x_train2[, j]))
    if(length(idx1) > 0) x_train2[idx1, j] <- mu
    
    idx2 <- which(is.na(x_new2[, j]))
    if(length(idx2) > 0) x_new2[idx2, j] <- mu
  }
  
  return(list(
    x_train = x_train2,
    x_new = x_new2,
    train_means = train_means
  ))
}

zscore_by_train <- function(x_train, x_new){
  x_train2 <- x_train
  x_new2   <- x_new
  
  train_means <- sapply(x_train2, function(z) mean(as.numeric(z), na.rm = TRUE))
  train_sds   <- sapply(x_train2, function(z) sd(as.numeric(z), na.rm = TRUE))
  
  train_sds[is.na(train_sds) | train_sds == 0] <- 1
  
  for(j in seq_len(ncol(x_train2))){
    x_train2[, j] <- (as.numeric(x_train2[, j]) - train_means[j]) / train_sds[j]
    x_new2[, j]   <- (as.numeric(x_new2[, j])   - train_means[j]) / train_sds[j]
  }
  
  return(list(
    x_train = as.data.frame(x_train2),
    x_new   = as.data.frame(x_new2),
    train_means = train_means,
    train_sds = train_sds
  ))
}

calc_clinical_metrics <- function(y_true, prob, thr){
  pred_class <- ifelse(prob >= thr, 1, 0)
  
  TP <- sum(pred_class == 1 & y_true == 1, na.rm = TRUE)
  TN <- sum(pred_class == 0 & y_true == 0, na.rm = TRUE)
  FP <- sum(pred_class == 1 & y_true == 0, na.rm = TRUE)
  FN <- sum(pred_class == 0 & y_true == 1, na.rm = TRUE)
  
  sens <- ifelse((TP + FN) == 0, NA, TP / (TP + FN))
  spec <- ifelse((TN + FP) == 0, NA, TN / (TN + FP))
  ppv  <- ifelse((TP + FP) == 0, NA, TP / (TP + FP))
  npv  <- ifelse((TN + FN) == 0, NA, TN / (TN + FN))
  acc  <- ifelse((TP + TN + FP + FN) == 0, NA, (TP + TN) / (TP + TN + FP + FN))
  
  precision <- ppv
  f1_score <- ifelse(
    is.na(precision) | is.na(sens) | (precision + sens) == 0,
    NA,
    2 * precision * sens / (precision + sens)
  )
  
  data.frame(
    TP = TP, TN = TN, FP = FP, FN = FN,
    Sensitivity = sens,
    Specificity = spec,
    PPV = ppv,
    Precision = precision,
    NPV = npv,
    Accuracy = acc,
    F1_score = f1_score
  )
}

# =======================
# 训练集预处理（参数只从训练集来）
# =======================
x_train_use <- x_train_raw
x_test_use  <- x_test_raw

if(do_impute_by_train_mean){
  imp_res <- impute_by_train_mean(x_train_use, x_test_use)
  x_train_use <- imp_res$x_train
  x_test_use  <- imp_res$x_new
  
  write.csv(
    data.frame(
      Feature = names(imp_res$train_means),
      Train_Mean = as.numeric(imp_res$train_means)
    ),
    "Train_mean_for_imputation.csv",
    row.names = FALSE
  )
}

# 去掉训练集方差为 0 的列
var0_cols <- sapply(x_train_use, function(z){
  z <- as.numeric(z)
  stats::var(z, na.rm = TRUE) == 0
})
var0_cols[is.na(var0_cols)] <- TRUE

if(any(var0_cols)){
  keep_cols <- names(var0_cols)[!var0_cols]
  x_train_use <- x_train_use[, keep_cols, drop = FALSE]
  x_test_use  <- x_test_use[,  keep_cols, drop = FALSE]
}

if(ncol(x_train_use) == 0) stop("过滤后没有可用特征列")

if(do_zscore_by_train){
  scale_res <- zscore_by_train(x_train_use, x_test_use)
  x_train_use <- scale_res$x_train
  x_test_use  <- scale_res$x_new
  
  write.csv(
    data.frame(
      Feature = names(scale_res$train_means),
      Train_Mean = as.numeric(scale_res$train_means),
      Train_SD   = as.numeric(scale_res$train_sds)
    ),
    "Train_zscore_parameters.csv",
    row.names = FALSE
  )
}

# =======================
# 建模：只用训练集
# =======================
train_model_df <- as.data.frame(x_train_use)
train_model_df$y <- y_train

logit <- glm(y ~ ., family = binomial(link = "logit"), data = train_model_df)

# 训练集预测概率
pred_train <- predict(logit, newdata = as.data.frame(x_train_use), type = "response")

# =======================
# 训练集 ROC + 全局阈值（只从训练集确定）
# =======================
roc_train <- pROC::roc(
  y_train,
  as.numeric(pred_train),
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

ci_train <- as.numeric(
  pROC::ci.auc(roc_train, method = "bootstrap", boot.n = boot_n_auc)
)

# 训练集 Youden 全局阈值
best_thr_info <- pROC::coords(
  roc_train,
  x = "best",
  best.method = "youden",
  ret = c("threshold", "sensitivity", "specificity"),
  transpose = FALSE
)

global_threshold <- as.numeric(best_thr_info["threshold"])

write.csv(
  data.frame(
    Threshold = global_threshold,
    Train_Sensitivity_at_best = as.numeric(best_thr_info["sensitivity"]),
    Train_Specificity_at_best = as.numeric(best_thr_info["specificity"])
  ),
  "Global_threshold_from_train.csv",
  row.names = FALSE
)

# 训练集临床指标
clin_train <- calc_clinical_metrics(
  y_true = y_train,
  prob   = as.numeric(pred_train),
  thr    = global_threshold
)

# 训练集预测输出
train_pred_out <- data.frame(
  Sample = rownames(x_train_use),
  Type = y_train,
  Pred_Prob = as.numeric(pred_train),
  Pred_Class = ifelse(as.numeric(pred_train) >= global_threshold, 1, 0),
  Threshold_Used = global_threshold,
  stringsAsFactors = FALSE
)
write.csv(train_pred_out, "Train_predictions.csv", row.names = FALSE)

# =======================
# 各外部数据集验证
# =======================
all_metrics_list <- list()

for(subset_name in names(test_groups)){
  idx <- test_groups[[subset_name]]
  
  x_test_subset <- x_test_use[idx, , drop = FALSE]
  y_test_subset <- y_test_all[idx]
  
  pred_test <- predict(logit, newdata = as.data.frame(x_test_subset), type = "response")
  
  roc_test <- pROC::roc(
    y_test_subset,
    as.numeric(pred_test),
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  ci_test <- as.numeric(
    pROC::ci.auc(roc_test, method = "bootstrap", boot.n = boot_n_auc)
  )
  
  clin_test <- calc_clinical_metrics(
    y_true = y_test_subset,
    prob   = as.numeric(pred_test),
    thr    = global_threshold
  )
  
  pred_out <- data.frame(
    Sample = rownames(x_test_subset),
    Type = y_test_subset,
    Pred_Prob = as.numeric(pred_test),
    Pred_Class = ifelse(as.numeric(pred_test) >= global_threshold, 1, 0),
    Threshold_Used = global_threshold,
    Dataset = subset_name,
    stringsAsFactors = FALSE
  )
  write.csv(pred_out, paste0(subset_name, "_predictions.csv"), row.names = FALSE)
  
  metrics_out <- data.frame(
    Cohort = subset_name,
    AUC = as.numeric(pROC::auc(roc_test)),
    CI_low = ci_test[1],
    CI_high = ci_test[3],
    Threshold = global_threshold,
    clin_test,
    stringsAsFactors = FALSE
  )
  
  all_metrics_list[[subset_name]] <- metrics_out
}

# =======================
# 汇总输出
# =======================
metrics_train <- data.frame(
  Cohort = "Train",
  AUC = as.numeric(pROC::auc(roc_train)),
  CI_low = ci_train[1],
  CI_high = ci_train[3],
  Threshold = global_threshold,
  clin_train,
  stringsAsFactors = FALSE
)

metrics_external <- do.call(rbind, all_metrics_list)
metrics_all <- rbind(metrics_train, metrics_external)

write.csv(metrics_all, "ROC.metrics.clinical.by_dataset.csv", row.names = FALSE)
print(metrics_all)

# =======================
# 汇总柱状图
# =======================
plot_metrics <- c("Sensitivity", "Specificity", "Precision", "Accuracy", "F1_score")

plot_df <- metrics_all[, c("Cohort", plot_metrics), drop = FALSE]

bar_long <- data.frame(
  Cohort = rep(plot_df$Cohort, each = length(plot_metrics)),
  Metric = rep(plot_metrics, times = nrow(plot_df)),
  Value  = as.vector(t(as.matrix(plot_df[, plot_metrics, drop = FALSE]))),
  stringsAsFactors = FALSE
)
bar_long <- bar_long[!is.na(bar_long$Value), , drop = FALSE]

pdf("Clinical_metrics_by_dataset.barplot.pdf", width = 12, height = 5)

ggplot(bar_long, aes(x = Metric, y = Value, fill = Cohort)) +
  geom_col(position = "dodge", width = 0.7) +
  coord_cartesian(ylim = c(0, 1)) +
  labs(
    title = paste0("Clinical metrics by dataset (global threshold from train = ",
                   round(global_threshold, 4), ")"),
    x = "",
    y = "Value"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

dev.off()

# =======================
# 训练集 + 外部验证 ROC 图
# =======================
pdf("ROC_all_datasets.pdf", width = 7, height = 7)

plot(roc_train, col = "#D55E00", lwd = 2, main = "ROC curves")
legend_text <- c(
  paste0("Train, AUC=", round(as.numeric(pROC::auc(roc_train)), 3))
)

legend_cols <- c("#D55E00")

for(i in seq_along(names(test_groups))){
  subset_name <- names(test_groups)[i]
  idx <- test_groups[[subset_name]]
  
  pred_test <- predict(logit,
                       newdata = as.data.frame(x_test_use[idx, , drop = FALSE]),
                       type = "response")
  
  roc_test <- pROC::roc(
    y_test_all[idx],
    as.numeric(pred_test),
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )
  
  this_col <- grDevices::rainbow(length(test_groups))[i]
  plot(roc_test, add = TRUE, col = this_col, lwd = 2)
  
  legend_text <- c(
    legend_text,
    paste0(subset_name, ", AUC=", round(as.numeric(pROC::auc(roc_test)), 3))
  )
  legend_cols <- c(legend_cols, this_col)
}

legend("bottomright", legend = legend_text, col = legend_cols, lwd = 2, cex = 0.9)
dev.off()