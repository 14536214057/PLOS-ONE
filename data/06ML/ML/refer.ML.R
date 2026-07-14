# =========================
# refer.ML.R
# =========================

RunML <- function(method, Train_set, Train_label, mode = "Model", classVar){
  method = gsub(" ", "", method)
  method_name  = gsub("(\\w+)\\[(.+)\\]", "\\1", method)
  method_param = gsub("(\\w+)\\[(.+)\\]", "\\2", method)
  
  method_param = switch(
    EXPR = method_name,
    "Enet"    = list("alpha" = as.numeric(gsub("alpha=", "", method_param))),
    "Stepglm" = list("direction" = method_param),
    NULL
  )
  
  message("Run ", method_name, " algorithm for ", mode, "; ",
          method_param, ";",
          " using ", ncol(Train_set), " Variables")
  
  args = list(
    "Train_set"   = Train_set,
    "Train_label" = Train_label,
    "mode"        = mode,
    "classVar"    = classVar
  )
  args = c(args, method_param)
  
  obj <- do.call(what = paste0("Run", method_name), args = args)
  
  if(mode == "Variable"){
    message(length(obj), " Variables retained;\n")
  }else{
    message("\n")
  }
  return(obj)
}

RunEnet <- function(Train_set, Train_label, mode, classVar, alpha){
  cv.fit = glmnet::cv.glmnet(
    x = Train_set,
    y = Train_label[[classVar]],
    family = "binomial",
    alpha = alpha,
    nfolds = 10
  )
  fit = glmnet::glmnet(
    x = Train_set,
    y = Train_label[[classVar]],
    family = "binomial",
    alpha = alpha,
    lambda = cv.fit$lambda.min
  )
  fit$subFeature = colnames(Train_set)
  if (mode == "Model")    return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLasso <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 1)
}

RunRidge <- function(Train_set, Train_label, mode, classVar){
  RunEnet(Train_set, Train_label, mode, classVar, alpha = 0)
}

RunStepglm <- function(Train_set, Train_label, mode, classVar, direction){
  fit <- stats::step(
    stats::glm(
      formula = Train_label[[classVar]] ~ .,
      family = "binomial",
      data = as.data.frame(Train_set)
    ),
    direction = direction,
    trace = 0
  )
  fit$subFeature = colnames(Train_set)
  if (mode == "Model")    return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunSVM <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit = e1071::svm(
    formula = eval(parse(text = paste(classVar, "~."))),
    data = data,
    probability = TRUE
  )
  fit$subFeature = colnames(Train_set)
  if (mode == "Model")    return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunLDA <- function(Train_set, Train_label, mode, classVar){
  data <- as.data.frame(Train_set)
  data[[classVar]] <- as.factor(Train_label[[classVar]])
  fit = caret::train(
    eval(parse(text = paste(classVar, "~."))),
    data = data,
    method = "lda",
    trControl = caret::trainControl(method = "cv")
  )
  fit$subFeature = colnames(Train_set)
  if (mode == "Model")    return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# ---- FIXED glmBoost (robust folds, avoid out-of-bounds) ----
RunglmBoost <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  
  fit <- mboost::glmboost(
    eval(parse(text = paste(classVar, "~."))),
    data = data,
    family = mboost::Binomial()
  )
  
  # ---- KEY FIX: weights must be length nrow(data) ----
  w <- rep(1, nrow(data))
  
  # stratified kfold to keep class ratio stable
  folds <- mboost::cv(w, type = "kfold", strata = data[[classVar]])
  
  cvm <- mboost::cvrisk(
    fit,
    papply = lapply,
    folds = folds
  )
  
  fit <- mboost::glmboost(
    eval(parse(text = paste(classVar, "~."))),
    data = data,
    family = mboost::Binomial(),
    control = mboost::boost_control(mstop = max(mboost::mstop(cvm), 40))
  )
  
  fit$subFeature = colnames(Train_set)
  if (mode == "Model") return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunplsRglm <- function(Train_set, Train_label, mode, classVar){
  nt_use <- min(10, ncol(Train_set))
  if (nt_use < 1) stop("plsRglm: Train_set has 0 variables.")
  
  suppressMessages(
    plsRglm::cv.plsRglm(
      formula = Train_label[[classVar]] ~ .,
      data = as.data.frame(Train_set),
      nt = nt_use,
      verbose = FALSE
    )
  )
  
  fit <- plsRglm::plsRglm(
    Train_label[[classVar]],
    as.data.frame(Train_set),
    modele = "pls-glm-logistic",
    nt = nt_use,
    verbose = FALSE,
    sparse = TRUE
  )
  
  fit$subFeature = colnames(Train_set)
  if (mode == "Model")    return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunRF <- function(Train_set, Train_label, mode, classVar){
  rf_nodesize = 5
  Train_label[[classVar]] <- as.factor(Train_label[[classVar]])
  
  fit <- randomForestSRC::rfsrc(
    formula = stats::formula(paste0(classVar, "~.")),
    data = cbind(Train_set, Train_label[classVar]),
    ntree = 1000,
    nodesize = rf_nodesize,
    importance = TRUE,
    proximity = TRUE,
    forest = TRUE
  )
  
  fit$subFeature = colnames(Train_set)
  if (mode == "Model")    return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunGBM <- function(Train_set, Train_label, mode, classVar){
  fit <- gbm::gbm(
    formula = Train_label[[classVar]] ~ .,
    data = as.data.frame(Train_set),
    distribution = "bernoulli",
    n.trees = 10000,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    cv.folds = 10,
    n.cores = 6
  )
  best <- which.min(fit$cv.error)
  
  fit <- gbm::gbm(
    formula = Train_label[[classVar]] ~ .,
    data = as.data.frame(Train_set),
    distribution = "bernoulli",
    n.trees = best,
    interaction.depth = 3,
    n.minobsinnode = 10,
    shrinkage = 0.001,
    n.cores = 8
  )
  
  fit$subFeature = colnames(Train_set)
  if (mode == "Model")    return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

# ---- FIXED XGBoost (force y to 0/1 for binary:logistic) ----
RunXGBoost <- function(Train_set, Train_label, mode, classVar){
  
  y_raw <- Train_label[[classVar]]
  
  # map to 0/1 safely
  if (is.factor(y_raw) || is.character(y_raw)) {
    y_fac <- as.factor(y_raw)
    if (nlevels(y_fac) != 2) stop("XGBoost: classVar must have exactly 2 classes.")
    y <- as.integer(y_fac) - 1
  } else {
    y_num <- as.numeric(y_raw)
    u <- sort(unique(y_num))
    if (length(u) != 2) stop("XGBoost: classVar must have exactly 2 unique numeric values.")
    y <- ifelse(y_num == u[1], 0, 1)
  }
  
  indexes = caret::createFolds(y, k = 5, list = TRUE)
  
  CV <- unlist(lapply(indexes, function(pt){
    dtrain = xgboost::xgb.DMatrix(
      data = Train_set[-pt, , drop = FALSE],
      label = y[-pt]
    )
    dtest = xgboost::xgb.DMatrix(
      data = Train_set[pt, , drop = FALSE],
      label = y[pt]
    )
    watchlist <- list(train = dtrain, test = dtest)
    
    bst <- xgboost::xgb.train(
      data = dtrain,
      max.depth = 2,
      eta = 1,
      nthread = 2,
      nrounds = 50,
      watchlist = watchlist,
      objective = "binary:logistic",
      eval_metric = "logloss",
      verbose = 0
    )
    which.min(bst$evaluation_log$test_logloss)
  }))
  
  nround <- as.numeric(names(which.max(table(CV))))
  
  fit <- xgboost::xgboost(
    data = Train_set,
    label = y,
    max.depth = 2,
    eta = 1,
    nthread = 2,
    nrounds = nround,
    objective = "binary:logistic",
    eval_metric = "logloss",
    verbose = 0
  )
  
  fit$subFeature = colnames(Train_set)
  if (mode == "Model")    return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

RunNaiveBayes <- function(Train_set, Train_label, mode, classVar){
  data <- cbind(Train_set, Train_label[classVar])
  data[[classVar]] <- as.factor(data[[classVar]])
  fit <- e1071::naiveBayes(eval(parse(text = paste(classVar, "~."))), data = data)
  fit$subFeature = colnames(Train_set)
  if (mode == "Model")    return(fit)
  if (mode == "Variable") return(ExtractVar(fit))
}

quiet <- function(..., messages=FALSE, cat=FALSE){
  if(!cat){
    sink(tempfile())
    on.exit(sink())
  }
  out <- if(messages) eval(...) else suppressMessages(eval(...))
  out
}

standarize.fun <- function(indata, centerFlag, scaleFlag){
  scale(indata, center=centerFlag, scale=scaleFlag)
}

scaleData <- function(data, cohort = NULL, centerFlags = NULL, scaleFlags = NULL){
  samplename = rownames(data)
  if (is.null(cohort)){
    data <- list(data); names(data) = "training"
  }else{
    data <- split(as.data.frame(data), cohort)
  }
  
  if (is.null(centerFlags)){
    centerFlags = FALSE; message("No centerFlags found, set as FALSE")
  }
  if (length(centerFlags) == 1){
    centerFlags = rep(centerFlags, length(data))
    message("set centerFlags for all cohort as ", unique(centerFlags))
  }
  if (is.null(names(centerFlags))){
    names(centerFlags) <- names(data)
    message("match centerFlags with cohort by order\n")
  }
  
  if (is.null(scaleFlags)){
    scaleFlags = FALSE; message("No scaleFlags found, set as FALSE")
  }
  if (length(scaleFlags) == 1){
    scaleFlags = rep(scaleFlags, length(data))
    message("set scaleFlags for all cohort as ", unique(scaleFlags))
  }
  if (is.null(names(scaleFlags))){
    names(scaleFlags) <- names(data)
    message("match scaleFlags with cohort by order\n")
  }
  
  centerFlags <- centerFlags[names(data)]
  scaleFlags  <- scaleFlags[names(data)]
  outdata <- mapply(standarize.fun, indata = data,
                    centerFlag = centerFlags, scaleFlag = scaleFlags,
                    SIMPLIFY = FALSE)
  outdata <- do.call(rbind, outdata)
  outdata <- outdata[samplename, ]
  return(outdata)
}

# ---- FIXED rfsrc feature extraction ----
ExtractVar <- function(fit){
  Feature <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"       = rownames(coef(fit))[which(coef(fit)[, 1] != 0)],
    "glm"          = names(coef(fit)),
    "svm.formula"  = fit$subFeature,
    "train"        = fit$coefnames,
    "glmboost"     = names(coef(fit)[abs(coef(fit)) > 0]),
    "plsRglmmodel" = rownames(fit$Coeffs)[fit$Coeffs != 0],
    
    "rfsrc" = {
      topv <- tryCatch(
        randomForestSRC::max.subtree(fit)$topvars,
        error = function(e) NULL
      )
      if (!is.null(topv) && length(topv) > 0) {
        topv
      } else {
        imp <- randomForestSRC::vimp(fit)$importance
        names(sort(imp, decreasing = TRUE))
      }
    },
    
    "gbm" = {
      s <- gbm::summary.gbm(fit, plotit = FALSE)
      rownames(s)[s$rel.inf > 0]
    },
    "xgb.Booster"  = fit$subFeature,
    "naiveBayes"   = fit$subFeature
  ))
  
  Feature <- setdiff(Feature, c("(Intercept)", "Intercept"))
  return(Feature)
}

CalPredictScore <- function(fit, new_data, type = "lp"){
  new_data <- new_data[, fit$subFeature, drop = FALSE]
  RS <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"       = predict(fit, type = "response", as.matrix(new_data)),
    "glm"          = predict(fit, type = "response", as.data.frame(new_data)),
    "svm.formula"  = predict(fit, as.data.frame(new_data), probability = TRUE),
    "train"        = predict(fit, new_data, type = "prob")[[2]],
    "glmboost"     = predict(fit, type = "response", as.data.frame(new_data)),
    "plsRglmmodel" = predict(fit, type = "response", as.data.frame(new_data)),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$predicted[, "1"],
    "gbm"          = predict(fit, type = "response", as.data.frame(new_data)),
    "xgb.Booster"  = predict(fit, as.matrix(new_data)),
    "naiveBayes"   = predict(object = fit, type = "raw", newdata = new_data)[, "1"]
  ))
  RS = as.numeric(as.vector(RS))
  names(RS) = rownames(new_data)
  return(RS)
}

PredictClass <- function(fit, new_data){
  new_data <- new_data[, fit$subFeature, drop = FALSE]
  label <- quiet(switch(
    EXPR = class(fit)[1],
    "lognet"       = predict(fit, type = "class", as.matrix(new_data)),
    "glm"          = ifelse(predict(fit, type = "response", as.data.frame(new_data)) > 0.5, "1", "0"),
    "svm.formula"  = predict(fit, as.data.frame(new_data), decision.values = TRUE),
    "train"        = predict(fit, new_data, type = "raw"),
    "glmboost"     = predict(fit, type = "class", as.data.frame(new_data)),
    "plsRglmmodel" = ifelse(predict(fit, type = "response", as.data.frame(new_data)) > 0.5, "1", "0"),
    "rfsrc"        = predict(fit, as.data.frame(new_data))$class,
    "gbm"          = ifelse(predict(fit, type = "response", as.data.frame(new_data)) > 0.5, "1", "0"),
    "xgb.Booster"  = ifelse(predict(fit, as.matrix(new_data)) > 0.5, "1", "0"),
    "naiveBayes"   = predict(object = fit, type = "class", newdata = new_data)
  ))
  label = as.character(as.vector(label))
  names(label) = rownames(new_data)
  return(label)
}

RunEval <- function(fit,
                    Test_set = NULL,
                    Test_label = NULL,
                    Train_set = NULL,
                    Train_label = NULL,
                    Train_name = NULL,
                    cohortVar = "Cohort",
                    classVar){
  
  if(!is.element(cohortVar, colnames(Test_label))) {
    stop(paste0("There is no [", cohortVar, "] indicator, please fill in one more column!"))
  }
  
  if((!is.null(Train_set)) & (!is.null(Train_label))) {
    new_data <- rbind.data.frame(
      Train_set[, fit$subFeature, drop = FALSE],
      Test_set[,  fit$subFeature, drop = FALSE]
    )
    
    if(!is.null(Train_name)) {
      Train_label$Cohort <- Train_name
    } else {
      Train_label$Cohort <- "Training"
    }
    colnames(Train_label)[ncol(Train_label)] <- cohortVar
    
    Test_label <- rbind.data.frame(
      Train_label[, c(cohortVar, classVar)],
      Test_label[,  c(cohortVar, classVar)]
    )
    Test_label[,1] <- factor(
      Test_label[,1],
      levels = c(unique(Train_label[,cohortVar]),
                 setdiff(unique(Test_label[,cohortVar]), unique(Train_label[,cohortVar])))
    )
  } else {
    new_data <- Test_set[, fit$subFeature, drop = FALSE]
  }
  
  RS <- suppressWarnings(CalPredictScore(fit = fit, new_data = new_data))
  
  Predict.out <- Test_label
  Predict.out$RS <- as.vector(RS)
  Predict.out <- split(x = Predict.out, f = Predict.out[,cohortVar])
  
  unlist(lapply(Predict.out, function(data){
    as.numeric(pROC::auc(suppressMessages(pROC::roc(data[[classVar]], data$RS))))
  }))
}

SimpleHeatmap <- function(Cindex_mat, avg_Cindex,
                          CohortCol, barCol,
                          cellwidth = 1, cellheight = 0.5,
                          cluster_columns, cluster_rows){
  
  col_ha = ComplexHeatmap::columnAnnotation(
    "Cohort" = colnames(Cindex_mat),
    col = list("Cohort" = CohortCol),
    show_annotation_name = FALSE
  )
  
  row_ha = ComplexHeatmap::rowAnnotation(
    bar = ComplexHeatmap::anno_barplot(
      avg_Cindex, bar_width = 0.8, border = FALSE,
      gp = grid::gpar(fill = barCol, col = NA),
      add_numbers = TRUE, numbers_offset = grid::unit(-10, "mm"),
      axis_param = list("labels_rot" = 0),
      numbers_gp = grid::gpar(fontsize = 9, col = "white"),
      width = grid::unit(3, "cm")
    ),
    show_annotation_name = FALSE
  )
  
  ComplexHeatmap::Heatmap(
    as.matrix(Cindex_mat), name = "AUC",
    right_annotation = row_ha,
    top_annotation = col_ha,
    col = c("#4195C1", "#FFFFFF", "#CB5746"),
    rect_gp = grid::gpar(col = "black", lwd = 1),
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = "left",
    width = grid::unit(cellwidth * ncol(Cindex_mat) + 2, "cm"),
    height = grid::unit(cellheight * nrow(Cindex_mat), "cm"),
    column_split = factor(colnames(Cindex_mat), levels = colnames(Cindex_mat)),
    column_title = NULL,
    cell_fun = function(j, i, x, y, w, h, col) {
      grid::grid.text(
        label = format(Cindex_mat[i, j], digits = 3, nsmall = 3),
        x, y, gp = grid::gpar(fontsize = 10)
      )
    }
  )
}

