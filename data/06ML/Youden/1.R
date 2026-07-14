library(pROC)
library(ggplot2)

# =======================
# 输入区：只改这里
# =======================
setwd("E:\\AAA_SCI\\IPF_A_SC\\10667_35145\\06ML\\Youden")

trainFile <- "data.train.txt"
testFile  <- "data.test.txt"

# threshold 计算方法：训练组 ROC 的 Youden 最优阈值
threshold_method <- "youden"

# 常规 baseline：不做 train/test 子集双向 Q3 靠近
do_q3_align_both <- FALSE

# 常规 baseline：不对 external/test 子集做异常值替换
do_outlier_replace_test <- FALSE

# 异常值阈值：超过 子集均值 ± outlier_sd_cutoff*SD 就替换
outlier_sd_cutoff <- 3

# bootstrap 次数
boot_n_auc <- 10000

# 从样本名提取数据集名称
# 当前 test 文件样本名格式如：
# GSE24206_GSM595407_Control
# GSE53845_GSM1302047_Control
extract_dataset_name <- function(x){
  sub("_.*$", "", x)
}

set.seed(123)

# =======================
# 读取 train/test（行=样本，列=特征+Type）
# =======================
train_df <- read.table(trainFile, header=TRUE, sep="\t", check.names=FALSE, row.names=1)
test_df  <- read.table(testFile,  header=TRUE, sep="\t", check.names=FALSE, row.names=1)

if(!("Type" %in% colnames(train_df))) stop("trainFile 缺少 Type 列")
if(!("Type" %in% colnames(test_df)))  stop("testFile 缺少 Type 列")

# =======================
# 分离 Type 和特征
# =======================
y_train <- as.numeric(as.character(train_df$Type))
y_test_all <- as.numeric(as.character(test_df$Type))

if(any(is.na(y_train))) stop("trainFile 的 Type 列无法正确转换为数值，请确认 Type 为 0/1")
if(any(is.na(y_test_all))) stop("testFile 的 Type 列无法正确转换为数值，请确认 Type 为 0/1")

x_train_raw <- train_df[, setdiff(colnames(train_df), "Type"), drop=FALSE]
x_test_raw  <- test_df[,  setdiff(colnames(test_df),  "Type"), drop=FALSE]

x_train_raw[] <- lapply(x_train_raw, function(z) as.numeric(as.character(z)))
x_test_raw[]  <- lapply(x_test_raw,  function(z) as.numeric(as.character(z)))

# 对齐共同特征
common_vars <- intersect(colnames(x_train_raw), colnames(x_test_raw))
if(length(common_vars) == 0) stop("train/test 没有共同特征列")

x_train_raw <- x_train_raw[, common_vars, drop=FALSE]
x_test_raw  <- x_test_raw[,  common_vars, drop=FALSE]

# 去掉全 NA 列
all_na_cols <- colnames(x_train_raw)[
  colSums(!is.na(x_train_raw)) == 0 | colSums(!is.na(x_test_raw)) == 0
]
if(length(all_na_cols) > 0){
  x_train_raw <- x_train_raw[, setdiff(colnames(x_train_raw), all_na_cols), drop=FALSE]
  x_test_raw  <- x_test_raw[,  setdiff(colnames(x_test_raw),  all_na_cols), drop=FALSE]
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

# 保存拆分信息
split_info <- data.frame(
  Sample = test_sample_names,
  Type = y_test_all,
  Dataset = test_dataset_names,
  stringsAsFactors = FALSE
)
write.csv(split_info, "Test_split_by_dataset.csv", row.names=FALSE)

# =======================
# 函数区
# =======================
replace_test_outliers_by_test_mean <- function(x_test, sd_cutoff=3){
  x_test_new <- x_test
  
  replace_n <- integer(ncol(x_test_new))
  names(replace_n) <- colnames(x_test_new)
  
  for(j in seq_len(ncol(x_test_new))){
    te <- as.numeric(x_test_new[, j])
    mu <- mean(te, na.rm=TRUE)
    s  <- sd(te, na.rm=TRUE)
    
    if(is.na(mu) || is.na(s) || s == 0) next
    
    lower <- mu - sd_cutoff * s
    upper <- mu + sd_cutoff * s
    
    idx <- which(!is.na(te) & (te < lower | te > upper))
    if(length(idx) > 0){
      x_test_new[idx, j] <- mu
      replace_n[j] <- length(idx)
    }
  }
  
  replace_info <- data.frame(
    Feature = names(replace_n),
    Replaced_N = as.integer(replace_n),
    stringsAsFactors = FALSE
  )
  
  return(list(
    x_test_new = x_test_new,
    replace_info = replace_info
  ))
}

q3_align_train_test_towards_center <- function(x_train, x_test){
  q3_train <- apply(x_train, 1, function(v) as.numeric(quantile(v, 0.75, na.rm=TRUE)))
  q3_test  <- apply(x_test,  1, function(v) as.numeric(quantile(v, 0.75, na.rm=TRUE)))
  
  ref_train <- median(q3_train, na.rm=TRUE)
  ref_test  <- median(q3_test,  na.rm=TRUE)
  ref_q3    <- mean(c(ref_train, ref_test), na.rm=TRUE)
  
  q3_train[q3_train == 0 | is.na(q3_train)] <- 1
  q3_test[q3_test == 0 | is.na(q3_test)] <- 1
  
  sf_train <- ref_q3 / q3_train
  sf_test  <- ref_q3 / q3_test
  
  x_train_aligned <- x_train * sf_train
  x_test_aligned  <- x_test  * sf_test
  
  return(list(
    x_train_aligned = x_train_aligned,
    x_test_aligned  = x_test_aligned,
    ref_train = ref_train,
    ref_test = ref_test,
    ref_q3 = ref_q3
  ))
}

impute_by_train_mean <- function(x_train, x_new){
  x_train2 <- x_train
  x_new2   <- x_new
  
  for(j in seq_len(ncol(x_train2))){
    mu <- mean(as.numeric(x_train2[, j]), na.rm=TRUE)
    if(is.na(mu)) mu <- 0
    
    idx1 <- which(is.na(x_train2[, j]))
    if(length(idx1) > 0) x_train2[idx1, j] <- mu
    
    idx2 <- which(is.na(x_new2[, j]))
    if(length(idx2) > 0) x_new2[idx2, j] <- mu
  }
  
  return(list(x_train = x_train2, x_new = x_new2))
}

get_best_threshold <- function(roc_obj, method="youden"){
  coords_res <- pROC::coords(
    roc_obj,
    x = "best",
    best.method = method,
    ret = c("threshold", "sensitivity", "specificity"),
    transpose = FALSE
  )
  
  if(is.data.frame(coords_res)){
    threshold_value <- as.numeric(coords_res$threshold[1])
    sens_value <- as.numeric(coords_res$sensitivity[1])
    spec_value <- as.numeric(coords_res$specificity[1])
  } else {
    threshold_value <- as.numeric(coords_res["threshold"])
    sens_value <- as.numeric(coords_res["sensitivity"])
    spec_value <- as.numeric(coords_res["specificity"])
  }
  
  return(list(
    threshold = threshold_value,
    sensitivity = sens_value,
    specificity = spec_value
  ))
}

calc_clinical_metrics <- function(y_true, prob, thr){
  pred_class <- ifelse(prob >= thr, 1, 0)
  
  TP <- sum(pred_class == 1 & y_true == 1, na.rm=TRUE)
  TN <- sum(pred_class == 0 & y_true == 0, na.rm=TRUE)
  FP <- sum(pred_class == 1 & y_true == 0, na.rm=TRUE)
  FN <- sum(pred_class == 0 & y_true == 1, na.rm=TRUE)
  
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

analyze_one_test_subset <- function(
    x_train_raw, y_train,
    x_test_subset_raw, y_test_subset,
    subset_name,
    do_q3_align_both = FALSE,
    do_outlier_replace_test = FALSE,
    outlier_sd_cutoff = 3,
    boot_n_auc = 10000,
    threshold_method = "youden"
){
  x_test_step <- x_test_subset_raw
  
  if(do_outlier_replace_test){
    outlier_res <- replace_test_outliers_by_test_mean(
      x_test = x_test_step,
      sd_cutoff = outlier_sd_cutoff
    )
    x_test_step <- outlier_res$x_test_new
    write.csv(
      outlier_res$replace_info,
      paste0(subset_name, "_outlier_replacement_summary.csv"),
      row.names = FALSE
    )
  }
  
  if(do_q3_align_both){
    q3_res <- q3_align_train_test_towards_center(x_train_raw, x_test_step)
    x_train_use <- q3_res$x_train_aligned
    x_test_use  <- q3_res$x_test_aligned
    
    q3_info <- data.frame(
      Cohort = subset_name,
      Train_Q3_Median = q3_res$ref_train,
      Test_Q3_Median  = q3_res$ref_test,
      Shared_Q3_Ref   = q3_res$ref_q3,
      stringsAsFactors = FALSE
    )
    write.csv(
      q3_info,
      paste0(subset_name, "_Q3_alignment_summary.csv"),
      row.names = FALSE
    )
  } else {
    x_train_use <- x_train_raw
    x_test_use  <- x_test_step
  }
  
  imp_res <- impute_by_train_mean(x_train_use, x_test_use)
  x_train_use <- imp_res$x_train
  x_test_use  <- imp_res$x_new
  
  var0_cols <- sapply(x_train_use, function(z){
    z <- as.numeric(z)
    stats::var(z, na.rm=TRUE) == 0
  })
  var0_cols[is.na(var0_cols)] <- TRUE
  
  if(any(var0_cols)){
    keep_cols <- names(var0_cols)[!var0_cols]
    x_train_use <- x_train_use[, keep_cols, drop=FALSE]
    x_test_use  <- x_test_use[,  keep_cols, drop=FALSE]
  }
  
  if(ncol(x_train_use) == 0) stop(paste0(subset_name, ": 过滤后没有可用特征列"))
  
  train_m <- as.data.frame(x_train_use)
  train_m$y <- y_train
  
  logit <- glm(y ~ ., family=binomial(link="logit"), data=train_m)
  
  pred_train <- predict(logit, newdata=as.data.frame(x_train_use), type="response")
  pred_test  <- predict(logit, newdata=as.data.frame(x_test_use),  type="response")
  
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
  
  # ===== 在训练组上计算最优 threshold =====
  thr_res <- get_best_threshold(roc_train, method = threshold_method)
  train_threshold <- thr_res$threshold
  
  threshold_info <- data.frame(
    Cohort = subset_name,
    Threshold_Method = threshold_method,
    Train_Threshold = train_threshold,
    Train_Sensitivity_at_Threshold = thr_res$sensitivity,
    Train_Specificity_at_Threshold = thr_res$specificity,
    stringsAsFactors = FALSE
  )
  write.csv(
    threshold_info,
    paste0(subset_name, "_train_threshold_summary.csv"),
    row.names = FALSE
  )
  
  clin_train <- calc_clinical_metrics(
    y_true = y_train,
    prob   = as.numeric(pred_train),
    thr    = train_threshold
  )
  
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
  
  # ===== 将训练组 threshold 应用到验证组 =====
  clin_test <- calc_clinical_metrics(
    y_true = y_test_subset,
    prob   = as.numeric(pred_test),
    thr    = train_threshold
  )
  
  metrics_out <- rbind(
    data.frame(
      Cohort = "Internal",
      Compared_with = subset_name,
      AUC = as.numeric(pROC::auc(roc_train)),
      CI_low = ci_train[1],
      CI_high = ci_train[3],
      Threshold = train_threshold,
      Threshold_Method = threshold_method,
      Threshold_Source = "Train",
      clin_train,
      stringsAsFactors = FALSE
    ),
    data.frame(
      Cohort = subset_name,
      Compared_with = subset_name,
      AUC = as.numeric(pROC::auc(roc_test)),
      CI_low = ci_test[1],
      CI_high = ci_test[3],
      Threshold = train_threshold,
      Threshold_Method = threshold_method,
      Threshold_Source = "Train",
      clin_test,
      stringsAsFactors = FALSE
    )
  )
  
  pred_out <- data.frame(
    Sample = rownames(x_test_subset_raw),
    Type = y_test_subset,
    Pred_Prob = as.numeric(pred_test),
    Pred_Class = ifelse(as.numeric(pred_test) >= train_threshold, 1, 0),
    Threshold_Used = train_threshold,
    Threshold_Source = "Train",
    stringsAsFactors = FALSE
  )
  
  write.csv(pred_out, paste0(subset_name, "_predictions.csv"), row.names = FALSE)
  
  return(list(
    metrics = metrics_out,
    pred = pred_out,
    threshold = train_threshold
  ))
}

# =======================
# 逐个数据集分析
# =======================
all_metrics_list <- list()
all_threshold_list <- list()
dataset_names <- names(test_groups)

for(i in seq_along(dataset_names)){
  subset_name <- dataset_names[i]
  idx <- test_groups[[i]]
  
  x_test_subset <- x_test_raw[idx, , drop=FALSE]
  y_test_subset <- y_test_all[idx]
  
  res_i <- analyze_one_test_subset(
    x_train_raw = x_train_raw,
    y_train = y_train,
    x_test_subset_raw = x_test_subset,
    y_test_subset = y_test_subset,
    subset_name = subset_name,
    do_q3_align_both = do_q3_align_both,
    do_outlier_replace_test = do_outlier_replace_test,
    outlier_sd_cutoff = outlier_sd_cutoff,
    boot_n_auc = boot_n_auc,
    threshold_method = threshold_method
  )
  
  all_metrics_list[[i]] <- res_i$metrics
  all_threshold_list[[i]] <- data.frame(
    Dataset = subset_name,
    Train_Threshold = res_i$threshold,
    Threshold_Method = threshold_method,
    stringsAsFactors = FALSE
  )
}

metrics_clinical_split <- do.call(rbind, all_metrics_list)
write.csv(metrics_clinical_split, "ROC.metrics.clinical.by_dataset.csv", row.names=FALSE)
print(metrics_clinical_split)

threshold_summary_all <- do.call(rbind, all_threshold_list)
write.csv(threshold_summary_all, "Train_thresholds_by_dataset.csv", row.names=FALSE)
print(threshold_summary_all)

# =======================
# 汇总柱状图
# =======================
plot_metrics <- c("Sensitivity", "Specificity", "Precision", "Accuracy", "F1_score")

plot_df <- metrics_clinical_split[, c("Cohort", "Compared_with", plot_metrics), drop=FALSE]

internal_first <- plot_df[plot_df$Cohort == "Internal", , drop=FALSE]
external_only  <- plot_df[plot_df$Cohort != "Internal", , drop=FALSE]

if(nrow(internal_first) > 0){
  plot_df2 <- rbind(internal_first[1, , drop=FALSE], external_only)
} else {
  plot_df2 <- external_only
}

bar_long <- data.frame(
  Cohort = rep(plot_df2$Cohort, each=length(plot_metrics)),
  Metric = rep(plot_metrics, times=nrow(plot_df2)),
  Value  = as.vector(t(as.matrix(plot_df2[, plot_metrics, drop=FALSE]))),
  stringsAsFactors = FALSE
)
bar_long <- bar_long[!is.na(bar_long$Value), , drop=FALSE]

pdf("Clinical_metrics_by_dataset.barplot.pdf", width=12, height=5)

ggplot(bar_long, aes(x=Metric, y=Value, fill=Cohort)) +
  geom_col(position="dodge", width=0.7) +
  coord_cartesian(ylim=c(0, 1)) +
  labs(
    title = paste0("Clinical metrics by dataset (train-derived threshold, method = ", threshold_method, ")"),
    x = "",
    y = "Value"
  ) +
  theme_bw() +
  theme(plot.title = element_text(hjust=0.5))

dev.off()