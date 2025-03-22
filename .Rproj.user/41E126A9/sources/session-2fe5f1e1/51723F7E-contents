






#'Functional_differential
#'
#' @param data Matrix of preprocessed single-cell metabolomics data
#' @param method "t_test", "ANOVA", "FC", "PLS-DA", "OPLS-DA", "RF_RFE", "Kruskal_Wallis", "svmrfeFeatureRanking"
#' @param typeFDR The method for false discovery rate (FDR) adjustment of p - values, with "BH" (Benjamini - Hochberg) as a common option.
#' @param algo The statistical test algorithm to be used, such as "t" for the t - test, to calculate the test statistic.
#' @param q The significance level threshold for the adjusted p - values (FDR), typically set to 0.05 to determine significant results.
#' @param ntree The number of trees to grow in the random forest model (default 500), affecting model stability and accuracy.
#' @param ntreeIterat The number of trees grown in each iteration of the variable - selection algorithm (default 300).
#' @param vars.drop.frac vars.drop.frac is the fraction (default 0.2) of variables dropped in each iteration of random forest variable selection, aiding in finding relevant variables.
#' @param cValue It is a parameter for the SVM (Support Vector Machine) in the SVM - RFE (Recursive Feature Elimination) process. Specifically, cValue (default is 50) represents the cost parameter in the SVM. A larger cValue makes the SVM try harder to classify all training samples correctly, potentially leading to overfitting, while a smaller value allows for more misclassifications in training to achieve better generalization.
#' @param Scaled It is a logical parameter indicating whether the input data should be scaled. When Scaled = TRUE, the input data will be scaled; when Scaled = FALSE (default), no scaling operation will be performed on the input data. Scaling can sometimes improve the performance of SVM by making features on a similar scale.
#' @param label_col A vector indicating the group labels (e.g., 0 and 1) for each sample in the data, used to distinguish different groups for the statistical test.
#'
#' @return A data frame or relevant result object depending on the selected differential analysis method. For example, if "t_test" is selected, it returns a data frame with columns "probeID", "Stat", "RawpValue", "AdjpValue".
#' @export
#'
#' @examples
#'\donttest{
#' }
#'
#'
differential <- function(data, label_col, method = c("t_test", "ANOVA", "FC", "PLS-DA", "OPLS-DA", "RF_RFE", "Kruskal_Wallis", "svmrfeFeatureRanking"),
                         typeFDR = "BH", algo = "t", q = 0.05, ntree = 500, ntreeIterat = 300, vars.drop.frac = 0.2, cValue = 50, Scaled = FALSE) {
  method <- match.arg(method)

  data2 <- apply(data[, -c(1:4)], 2, as.numeric)
  data4 <- data[, 1:4]

  data2[data2 == 0] <- NA

  labels <- as.numeric(as.factor(data[, label_col]))


  require(multtest)
  require(mixOmics)
  require(e1071)
  library(varSelRF)
  library(ropls)


  T_test <- function(data, labels, typeFDR = "BH", algo = "t", q = 0.05) {
    list.probe <- rownames(data)
    if (is.null(list.probe)) list.probe <- 1:nrow(data)

    message("Launch ", algo, " test")
    test <- mt.teststat(data, labels, test = algo)

    s1 <- apply(data[, which(labels == 0)], 1, function(x) length(which(!is.na(x))))
    s2 <- apply(data[, which(labels == 1)], 1, function(x) length(which(!is.na(x))))

    w1 <- apply(data[, which(labels == 0)], 1, var, na.rm = TRUE) / s1
    w2 <- apply(data[, which(labels == 1)], 1, var, na.rm = TRUE) / s2

    ddl <- (w1 + w2)^2 / ((w1^2 / (s1 - 1)) + (w2^2 / (s2 - 1)))
    pval.test <- 2 * (pt(-abs(test), ddl))

    padj <- p.adjust(pval.test, method = typeFDR)
    out <- data.frame(list.probe, test, pval.test, padj)
    colnames(out) <- c("probeID", "Stat", "RawpValue", "AdjpValue")

    return(out)
  }

  one_ANOVA <- function(data, labels) {
    data <- as.data.frame(t(data))
    ANOVA.res <- NULL

    for(i in 1:dim(data)[2]){
      data_la <- cbind(data, labels)
      ANOVA.test <- aov(data[,i]~labels,data=data_la)
      ANOVA.s <- summary(ANOVA.test)[[1]]["Pr(>F)"][[1]][1]
      p.res <- c(colnames(data)[i], ANOVA.s)
      ANOVA.res <- rbind(ANOVA.res, p.res)
    }

    ANOVA.res <- as.data.frame(ANOVA.res)
    colnames(ANOVA.res) <- c("feature", "p.value")
    row.names(ANOVA.res) <- NULL

    return(ANOVA.res)
  }

  FC_test <- function(data, labels, method = mean, unlog = FALSE) {
    if (missing(data))
      stop("** FAILURE : 'data' is missing **")
    if (missing(labels))
      stop("** FAILURE : 'labels' is missing **")

    data1.ave <- apply(data[, which(labels == 1)], 1, method, na.rm = TRUE)
    data2.ave <- apply(data[, which(labels == 2)], 1, method, na.rm = TRUE)
    fold.change = data2.ave - data1.ave

    result <- data.frame(Variable = names(fold.change), ExpreInControl = data1.ave,
                         ExpreInCase = data2.ave, FoldChange = round(fold.change, 5))
    rownames(result) <- NULL
    return(result)
  }

  PLSDA_test <- function(mat, label) {

    X <- t(mat)
    Y <- as.factor(label)


    plsda.res <- plsda(X, Y, ncomp = 2)


    VIP <- mixOmics::vip(plsda.res)
    tab <- as.data.frame(VIP[order(VIP[, ncol(VIP)], decreasing = TRUE), ncol(VIP)])
    colnames(tab) <- "VIP"


    result <- list(tab = tab, sup1 = rownames(tab)[which(tab$VIP > 1)])
    class(result) <- "PLSDA.VIP"
    return(result)
  }


  OPLSDA_test <- function(mat, lab) {


    X <- t(mat)
    Y <- as.numeric(as.factor(lab))

    oplsda <- opls(X, Y, predI = 1, orthoI = 2, fig.pdfC = "none")

    res <- oplsda@vipVn
    cpds <- data.frame(CompoundName = names(res), VIP = res)
    return(cpds)
  }

  RF_test <- function(mat, label, ntree = 500, ntreeIterat = 300,
                      vars.drop.frac = 0.2) {
    set.seed(3)
    x <- t(mat)
    cl <- factor(label)
    rf.vs1 <- varSelRF(x, cl, ntree = ntree, ntreeIterat = ntreeIterat, vars.drop.frac = vars.drop.frac)

    RF.test <- rf.vs1$initialImportances

    nu_l <- match(rf.vs1$selected.vars, row.names(RF.test))

    num <- 1:dim(RF.test)[1]

    num[nu_l] <- "TRUE"
    num[-nu_l] <- "FALSE"

    RF.result <- cbind(row.names(RF.test), RF.test, as.data.frame(num))

    colnames(RF.result)[1] <- "Feature"
    colnames(RF.result)[3] <- "Significance"

    return(RF.result)
  }

  Kruskal_Wallis <- function(data, labels) {
    if (missing(data))
      stop("** FAILURE : 'data' is missing **")
    if (missing(labels))
      stop("** FAILURE : 'labels' is missing **")

    kru.res <- NULL

    for(i in 1:dim(data)[2]){
      data_la <- cbind(data, labels)
      kru <- kruskal.test(data[,i]~labels,data=data_la)

      p.res <- c(colnames(data)[i], kru$p.value)

      kru.res <- rbind(kru.res, p.res)
    }

    kru.res <- as.data.frame(kru.res)
    colnames(kru.res) <- c("feature", "p.value")
    row.names(kru.res) <- NULL

    return(kru.res)
  }

  rsvmRFE.wrap <- function(test.fold, X, ...) {
    train.data = X[-test.fold, ]
    test.data  = X[test.fold, ]

    features.ranked = svmRFE(train.data, ...)

    return(list(feature.ids=features.ranked, train.data.ids=row.names(train.data), test.data.ids=row.names(test.data)))
  }

  svmRFE <- function(X, k=1, halve.above=5000) {
    n = ncol(X) - 1

    cat('Scaling data...')
    X[, -1] = scale(X[, -1])
    cat('Done!\n')
    flush.console()

    pb = txtProgressBar(1, n, 1, style=3)

    i.surviving = 1:n
    i.ranked    = n
    ranked.list = vector(length=n)

    while(length(i.surviving) > 0) {
      if(k > 1) {
        folds = rep(1:k, len=nrow(X))[sample(nrow(X))]
        folds = lapply(1:k, function(x) which(folds == x))

        w = lapply(folds, getWeights, X[, c(1, 1+i.surviving)])
        w = do.call(rbind, w)

        w = t(apply(w, 1, function(x) x / sqrt(sum(x^2))))

        v    = w * w
        vbar = apply(v, 2, mean)
        vsd  = apply(v, 2, sd)
        c    = vbar / vsd
      } else {
        w = getWeights(NULL, X[, c(1, 1+i.surviving)])
        c = w * w
      }

      ranking = sort(c, index.return=T)$ix
      if(length(i.surviving) == 1) {
        ranking = 1
      }

      if(length(i.surviving) > halve.above) {
        nfeat = length(i.surviving)
        ncut  = round(nfeat / 2)
        n     = nfeat - ncut

        cat('Features halved from', nfeat, 'to', n, '\n')
        flush.console()

        pb = txtProgressBar(1, n, 1, style=3)

      } else ncut = 1

      ranked.list[i.ranked:(i.ranked-ncut+1)] = i.surviving[ranking[1:ncut]]
      i.ranked    = i.ranked - ncut
      i.surviving = i.surviving[-ranking[1:ncut]]

      setTxtProgressBar(pb, n-length(i.surviving))
      flush.console()
    }

    close(pb)

    return (ranked.list)
  }

  getWeights <- function(test.fold, X) {
    train.data = X
    if(!is.null(test.fold)) train.data = X[-test.fold, ]

    svmModel = svm(train.data[, -1], train.data[, 1], cost=10, cachesize=500,
                   scale=F, type="C-classification", kernel="linear")

    t(svmModel$coefs) %*% svmModel$SV
  }

  WriteFeatures <- function(results, input, save=T, file='features_ranked.txt') {
    featureID = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$ix
    avg.rank  = sort(apply(sapply(results, function(x) sort(x$feature, index.return=T)$ix), 1, mean), index=T)$x
    feature.name = colnames(input[, -1])[featureID]
    features.ranked = data.frame(FeatureName=feature.name, FeatureID=featureID, AvgRank=avg.rank)
    if(save==T) {
      write.table(features.ranked, file=file, quote=F, row.names=F)
    } else {
      return(features.ranked)
    }
  }

  svmrfeFeatureRanking <- function(data, labels, cValue = 50, Scaled = FALSE){
    input <- cbind(labels, data)

    results <- lapply(1:5, svmRFE.wrap, input, k=10, halve.above=100)

    top.features <- WriteFeatures(results, input, save=F)

    return (top.features)
  }


  if (method == "t_test") {

    labels <- as.numeric(as.factor(labels))
    labels <- ifelse(labels == 2, 1, 0)
    T_res <- T_test(t(data2), labels)
    T_res1 <- T_res[, c("probeID", "RawpValue")]
    marker.table2 <- T_res1[order(T_res1$RawpValue, decreasing = FALSE), ]
    return(marker.table2)

  } else if (method == "ANOVA") {

    marker.table1 <- one_ANOVA(t(data2), labels)
    marker.table2 <- marker.table1[order(as.numeric(as.vector(marker.table1[,2]))), ]
    return(marker.table2)

  } else if (method == "FC") {

    FCres <- FC_test(t(data2), as.numeric(as.factor(labels)))
    marker_so <- FCres[, c("Variable", "FoldChange")]
    marker.table2 <- marker_so[order(marker_so$FoldChange, decreasing = TRUE), ]
    return(marker.table2)

  } else if (method == "PLS-DA") {

    PLSDA_res <- PLSDA_test(t(data2), as.numeric(labels))
    PLSDA_res$tab$name <- rownames(PLSDA_res$tab)
    PLSDA_res$tab <- PLSDA_res$tab[order(PLSDA_res$tab$VIP, decreasing = TRUE), ]
    marker.table1 <- PLSDA_res$tab
    marker.table2 <- marker.table1[, c(2, 1)]
    return(marker.table2)

  } else if (method == "OPLS-DA") {

    marker.table <- OPLSDA_test(t(data2), labels)
    marker.table2 <- marker.table[order(marker.table$VIP, decreasing = TRUE), ]
    return(marker.table2)

  } else if (method == "RF_RFE") {

    RF_test_res <- RF_test(t(data2), labels)
    RF_test_res_so <- RF_test_res[, c(1, 2)]
    marker.table2 <- RF_test_res_so[order(-RF_test_res_so$MeanDecreaseAccuracy), ]
    return(marker.table2)

  } else if (method == "Kruskal_Wallis") {

    return(Kruskal_Wallis(data2, labels))

  } else if (method == "svmrfeFeatureRanking") {

    return(svmrfeFeatureRanking(t(data2), labels, cValue, Scaled))

  }
}
