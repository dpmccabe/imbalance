sapply(unlist(strsplit("MASS MBESS dplyr plyr ggplot2 gridExtra DMwR doMC RANN unbalanced ROSE xtable klaR matrixStats randomForest", " ")), function(pkg) {
  if (!is.element(pkg, installed.packages()[, 1])) install.packages(pkg, dep = T)
  library(pkg, character.only = T, quietly = T)
})

rm(list = ls())
registerDoMC(cores = 16)
setMKLthreads(16)

######################################################################
# adasyn
######################################################################
adasyn <- function(form, d, k = 5, beta = 1, do_parallel = T) {
  outcome_colname <- all.vars(form)[1]
  outcome_col <- which(colnames(d) == outcome_colname)
  pred_colnames <- all.vars(form)[-1]
  npred <- length(pred_colnames)
  minor_clname <- levels(d[, outcome_colname])[which.min(table(d[, outcome_colname]))]
  major_clname <- setdiff(levels(d[, outcome_colname]), minor_clname)

  d2 <- cbind(minor = d[, outcome_colname] == minor_clname, d[, pred_colnames], synthetic = F)

  minor_needed <- (sum(!d2$minor) - sum(d2$minor)) * beta
  minor_ind <- which(d2$minor)
  new_xs <- data.frame()

  gammas <- unname(aaply(1:length(minor_ind), 1, function(i) {
    x <- d2[minor_ind[i], pred_colnames]
    x_nn_inds <- nn2(d2[, pred_colnames], x, k = k + 1)$nn.idx[, -1]
    gam <- sum(!d2[x_nn_inds, "minor"]) / k
    return(gam)
  }, .parallel = do_parallel))

  norm_gammas <- round(minor_needed * gammas / sum(gammas))
  to_synth_ind <- minor_ind[norm_gammas > 0]
  to_synth_n <- norm_gammas[norm_gammas > 0]

  if (length(to_synth_ind) > 0) {
    new_xs <- adply(1:length(to_synth_ind), 1, function(j) {
      synth_ref <- d2[to_synth_ind[j], pred_colnames]
      nn_inds <- nn2(d2[minor_ind, pred_colnames], synth_ref, k = k + 1)$nn.idx[, -1]

      new_xs_j <- adply(1:to_synth_n[j], 1, function(h) {
        return(synth_ref + (d2[minor_ind[sample(nn_inds, 1)], pred_colnames] - synth_ref) * runif(npred))
      })[, -1]

      return(cbind(new_xs_j, synthetic = T))
    }, .parallel = do_parallel)[, -1]

    d2 <- rbind(d2, data.frame(minor = T, new_xs))
  }

  d2[, outcome_colname] <- factor(ifelse(d2$minor, minor_clname, major_clname))

  return(subset(d2, select = -c(minor, synthetic)))
}

######################################################################
# bSMOTE
######################################################################
bSMOTE <- function(form, d, k = 5, s = 5, do_parallel = T) {
  outcome_colname <- all.vars(form)[1]
  outcome_col <- which(colnames(d) == outcome_colname)
  pred_colnames <- all.vars(form)[-1]
  npred <- length(pred_colnames)
  minor_clname <- levels(d[, outcome_colname])[which.min(table(d[, outcome_colname]))]
  major_clname <- setdiff(levels(d[, outcome_colname]), minor_clname)

  d2 <- cbind(minor = d[, outcome_colname] == minor_clname, d[, pred_colnames], synthetic = F)

  minor_ind <- which(d2$minor)
  new_xs <- data.frame()

  danger <- unname(aaply(1:length(minor_ind), 1, function(i) {
    x <- d2[minor_ind[i], pred_colnames]
    x_nn_inds <- nn2(d2[, pred_colnames], x, k = k + 1)$nn.idx[, -1]
    nn_major <- sum(!d2[x_nn_inds, "minor"])
    return(nn_major >= k / 2 & nn_major < k)
  }, .parallel = do_parallel))

  to_synth_ind <- which(danger)

  if (length(to_synth_ind) > 0) {
    new_xs <- adply(1:length(to_synth_ind), 1, function(j) {
      synth_ref <- d2[to_synth_ind[j], pred_colnames]
      nn_inds <- nn2(d2[minor_ind, pred_colnames], synth_ref, k = k + 1)$nn.idx[, -1]
      nn_inds_smp <- sample(nn_inds, s)

      new_xs_j <- adply(1:s, 1, function(h) {
        return(synth_ref + (d2[minor_ind[nn_inds_smp[h]], pred_colnames] - synth_ref) * runif(s))
      })[, -1]

      return(cbind(new_xs_j, synthetic = T))
    }, .parallel = do_parallel)[, -1]

    d2 <- rbind(d2, data.frame(minor = T, new_xs))
  }

  d2[, outcome_colname] <- factor(ifelse(d2$minor, minor_clname, major_clname))

  return(subset(d2, select = -c(minor, synthetic)))
}

######################################################################
# example data
######################################################################
set.seed(20150613)
ex_init1 <- data.frame(mvrnorm(1600, c(10, 10), cor2cov(matrix(c(1, 0.3, 0.3, 1), 2, 2), c(4, 4))))
ex_init2 <- data.frame(mvrnorm(400, c(15, 15), cor2cov(matrix(c(1, -0.7, -0.7, 1), 2, 2), c(1.2, 1.2))))
ex_init3 <- data.frame(mvrnorm(1200, c(5, 6), cor2cov(matrix(c(1, 0.1, 0.1, 1), 2, 2), c(2, 2))))
ex_init = rbind(ex_init1, ex_init2, ex_init3)
colnames(ex_init) = c("x1", "x2")

ex_minor_ref1_x <- mvrnorm(100, c(7, 12), cor2cov(matrix(c(1, 0.5, 0.5, 1), 2, 2), c(1.5, 0.75)))
ex_minor_ref2_x <- mvrnorm(50, c(13, 6), cor2cov(matrix(c(1, -0.2, -0.2, 1), 2, 2), c(0.5, 1)))

ex_init <- filter(ex_init, x1 >= 0, x2 >= 0, x1 <= 20, x2 <= 20)

ex_minor_nn <- nn2(ex_init, rbind(ex_minor_ref1_x, ex_minor_ref2_x), k = 5)$nn.idx
ex_major_switch_i <- apply(ex_minor_nn, 1, function(r) {
  sample(r, 2, prob = 5:1)
})
ex_major_switch_i2 <- sample.int(nrow(ex_init), 15, prob = 1 / (ex_init$x1 / ex_init$x2))
ex_major_switch_i <- unique(c(ex_major_switch_i, ex_major_switch_i2))

ex_minor <- data.frame(ex_init[ex_major_switch_i,], type = "A")
ex_major <- data.frame(ex_init[-ex_major_switch_i,], type = "B")
ex <- rbind(ex_minor, ex_major)
ex$id <- 1:nrow(ex)
ex$type <- relevel(ex$type, "A")

summary(ex)
table(ex$type) / sum(table(ex$type))
ggplot(ex, aes(x1, x2)) + geom_point(size = 1.2, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.9))

######################################################################
# bad classifiers
######################################################################
bad <- data.frame(pred = rep("B", nrow(ex)), obs = ex$type)
(bad_cm <- confusionMatrix(bad$pred, bad$obs, positive = "A"))
t(bad_cm$table)
roc.curve(bad$obs, bad$pred)

bad2 <- data.frame(pred = c(rep("A", 50), rep("B", nrow(ex) - 60), rep("A", 10)), obs = ex$type)
(bad2_cm <- confusionMatrix(bad2$pred, bad2$obs, positive = "A"))
t(bad2_cm$table)
roc.curve(bad2$obs, bad2$pred)

######################################################################
# more plots
######################################################################
grid.arrange(ggplot(filter(ex, type == "B"), aes(x1, x2)) + geom_point(size = 1.1, aes(color = type), shape = 17) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("blue")) + theme(legend.position = c(0.9, 0.9)), ggplot(filter(ex, type == "A"), aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red")) + theme(legend.position = c(0.9, 0.9)), ncol = 2)

ex_with_safety <- ddply(ex, .(id), function(r) {
  if (r$type == "B") {
    safety <- NA
  } else {
    nn <- nn2(ex[, 1:2], r, k = 5)$nn.idx[, -1]
    n_minor_nn <- sum(ex[nn, "type"] == "A")
    safety <- c("noise", "borderline", "borderline", "safe")[n_minor_nn + 1]
  }

  return(cbind(r, safety))
}, .parallel = T)

safe <- which(ex_with_safety$safety == "safe")
borderline <- which(ex_with_safety$safety == "borderline")
noise <- which(ex_with_safety$safety == "noise")

p1 <- ggplot(ex_with_safety, aes(x1, x2)) + geom_point(size = 1, aes(color = type, shape = type), alpha = 0.25) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.9)) + annotate("point", ex_with_safety[safe, "x1"], ex_with_safety[safe, "x2"], size = 1, color = "red") + theme(legend.position = "none")

p2 <- ggplot(ex_with_safety, aes(x1, x2)) + geom_point(size = 1, aes(color = type, shape = type), alpha = 0.25) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.9)) + annotate("point", ex_with_safety[borderline, "x1"], ex_with_safety[borderline, "x2"], size = 1, color = "red") + theme(legend.position = "none")

p3 <- ggplot(ex_with_safety, aes(x1, x2)) + geom_point(size = 1, aes(color = type, shape = type), alpha = 0.25) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.9)) + annotate("point", ex_with_safety[noise, "x1"], ex_with_safety[noise, "x2"], size = 1, color = "red") + theme(legend.position = "none")

grid.arrange(p1, p2, p3, ncol = 1)

######################################################################
# example resampling
######################################################################
p_ex_none <- ggplot(ex, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = "none")
p_ex_none_alpha <- ggplot(ex, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type), alpha = 0.25) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = "none")

# under
ex_ubunder <- ubUnder(ex[, 1:2], ifelse(ex$type == "A", 1, 0))
str(ex_ubunder)
ex_under <- cbind(ex_ubunder$X, type = ex_ubunder$Y)
ex_under$type <- as.factor(ifelse(ex_under$type == 1, "A", "B"))
p_ex_under <- ggplot(ex_under, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_under, ncol = 2)

# under (20%)
ex_ubunder2 <- ubUnder(ex[, 1:2], ifelse(ex$type == "A", 1, 0), perc = 20, method = "percPos")
str(ex_ubunder2)
ex_under2 <- cbind(ex_ubunder2$X, type = ex_ubunder2$Y)
ex_under2$type <- as.factor(ifelse(ex_under2$type == 1, "A", "B"))
p_ex_under2 <- ggplot(ex_under2, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_under2, ncol = 2)

# over
ex_ubover <- ubOver(ex[, 1:2], ifelse(ex$type == "A", 1, 0))
str(ex_ubover)
ex_over <- cbind(ex_ubover$X, type = ex_ubover$Y)
ex_over$type <- as.factor(ifelse(ex_over$type == 1, "A", "B"))
p_ex_over <- ggplot(ex_over, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type), alpha = 0.25) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none_alpha, p_ex_over, ncol = 2)

# SMOTE
ex_ubsmote <- ubSMOTE(ex[, 1:2], as.factor(ifelse(ex$type == "A", 1, 0)), perc.over = 200, perc.under = 1000)
str(ex_ubsmote)
ex_smote <- cbind(ex_ubsmote$X, type = ex_ubsmote$Y)
ex_smote$type <- as.factor(ifelse(ex_smote$type == 1, "A", "B"))
p_ex_smote <- ggplot(ex_smote, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_smote, ncol = 2)

# ADASYN
ex_adasyn <- adasyn(type ~ x1 + x2, ex, beta = 0.1)
summary(ex_adasyn)
p_ex_adasyn <- ggplot(ex_adasyn, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_adasyn, ncol = 2)

# ROSE
ex_ubrose <- ROSE(type ~ x1 + x2, ex, p = 0.2, hmult.mino = 0.2)
summary(ex_ubrose)
ex_rose <- ex_ubrose$data
ex_rose$type <- relevel(ex_rose$type, "A")
p_ex_rose <- ggplot(ex_rose, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type), na.rm = T) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_rose, ncol = 2)

# Borderline-SMOTE
ex_bsmote <- bSMOTE(type ~ x1 + x2, ex, k = 5, s = 5)
summary(ex_bsmote)
p_ex_bsmote <- ggplot(ex_bsmote, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_bsmote, ncol = 2)

# Tomek
ex_ubtomek <- ubTomek(ex[, 1:2], ifelse(ex$type == "A", 1, 0))
str(ex_ubtomek)
ex_tomek <- cbind(ex_ubtomek$X, type = ex_ubtomek$Y)
ex_tomek$type <- as.factor(ifelse(ex_tomek$type == 1, "A", "B"))
p_ex_tomek = ggplot(ex_tomek, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_tomek, ncol = 2)

# SMOTE+Tomek
ex_ubsmotetomek <- ubTomek(ex_smote[, 1:2], ifelse(ex_smote$type == "A", 1, 0))
str(ex_ubsmotetomek)
ex_smotetomek <- cbind(ex_ubsmotetomek$X, type = ex_ubsmotetomek$Y)
ex_smotetomek$type <- as.factor(ifelse(ex_smotetomek$type == 1, "A", "B"))
p_ex_smotetomek <- ggplot(ex_smotetomek, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_smotetomek, ncol = 2)

# Tomek+SMOTE
ex_ubtomeksmote <- ubSMOTE(ex_tomek[, 1:2], as.factor(ifelse(ex_tomek$type == "A", 1, 0)), perc.over = 200, perc.under = 1000)
str(ex_ubtomeksmote)
ex_tomeksmote <- cbind(ex_ubtomeksmote$X, type = ex_ubtomeksmote$Y)
ex_tomeksmote$type <- as.factor(ifelse(ex_tomeksmote$type == "1", "A", "B"))
p_ex_tomeksmote <- ggplot(ex_tomeksmote, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_tomeksmote, ncol = 2)

# NCL
ex_ubncl <- ubNCL(ex[, 1:2], ifelse(ex$type == "A", 1, 0))
str(ex_ubncl)
ex_ncl <- cbind(ex_ubncl$X, type = ex_ubncl$Y)
ex_ncl$type <- as.factor(ifelse(ex_ncl$type == "1", "A", "B"))
p_ex_ncl <- ggplot(ex_ncl, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_ncl, ncol = 2)

# bSMOTE+NCL
ex_ubbsmotencl <- ubNCL(ex_bsmote[, 1:2], ifelse(ex_bsmote$type == "A", 1, 0))
str(ex_ubbsmotencl)
ex_bsmotencl <- cbind(ex_ubbsmotencl$X, type = ex_ubbsmotencl$Y)
ex_bsmotencl$type <- as.factor(ifelse(ex_bsmotencl$type == "1", "A", "B"))
p_ex_bsmotencl <- ggplot(ex_bsmotencl, aes(x1, x2)) + geom_point(size = 1.1, aes(color = type, shape = type)) + xlim(c(0, 20)) + ylim(c(0, 20)) + scale_colour_manual(values = c("red", "blue")) + theme(legend.position = c(0.9, 0.85))
grid.arrange(p_ex_none, p_ex_bsmotencl, ncol = 2)

######################################################################
# load mediastinal
######################################################################
med_ids <- unname(unlist(c(read.csv("data/med_ids.csv", header = F, stringsAsFactors = F))))
med <- read.delim("data/Mediastinal.res", header = F, skip = 3, stringsAsFactors = F)
med <- med[, c(2, seq(1, ncol(med), by = 2))]
str(med)

medt <- as.data.frame(t(med[, -c(1, 2)]))
gene_names <- data.frame(v = colnames(medt), gene = med[, 1])

medt$id <- med_ids

med_classes <- read.delim("data/Mediastinal.txt")
med_classes$group <- sample(rep(1:16, length.out = nrow(med_classes)))
med2 <- merge(med_classes, medt, by.x = "Sample.ID", by.y = "id")
colnames(med2)[1:2] <- c("id", "type")
med2$type <- relevel(med2$type, ref = "MLBCL")
str(med2)
table(med2$type)

rm(list = c("med", "medt"))

mform <- formula(type ~ . - id - group)
med_preds <- colnames(med2)[-c(1, 2, 3)]

######################################################################
# compute MAD, SNR
######################################################################
mads <- apply(med2[, med_preds], 2, function(g) mean(abs(g - median(g))))
top_genes <- names(tail(sort(mads), 15000))

get_snr_genes <- function(d, cols, m) {
  mu1 <- colMeans(subset(d, type == "MLBCL", select = cols))
  mu0 <- colMeans(subset(d, type == "DLBCL", select = cols))
  sigma1 <- colSds(as.matrix(subset(d, type == "MLBCL", select = cols)))
  sigma0 <- colSds(as.matrix(subset(d, type == "DLBCL", select = cols)))

  snr <- (mu1 - mu0) / (sigma1 + sigma0)
  return(c(names(tail(sort(snr), m / 2)), names(tail(sort(snr, decreasing = T), m / 2))))
}

m <- 100

top_genes_loo <- dlply(med2, .(id), function(r) {
  med_train <- subset(med2, id != r$id)
  train_genes <- get_snr_genes(med_train, top_genes, m)
  return(train_genes)
}, .parallel = T)

######################################################################
# resampling comparison
######################################################################
resamp_methods <- expand.grid(algo = "rf", do_mad = F, resamp_when = c("after_vs" , "before_vs"), which = c("None", "Over2", "Over4", "Over0", "Under25", "Under50", "SMOTE_2_2", "SMOTE_3_1", "SMOTE_1_3", "Tomek", "NCL3", "NCL4", "ROSE", "ADASYN", "bSMOTE", "SMOTE_2_2_Tomek", "Tomek_SMOTE_2_2", "SMOTE_2_2_NCL3", "NCL3_SMOTE_2_2"))
resamp_methods <- resamp_methods[-2,]
resamp_methods <- rbind(expand.grid(algo = "nb", do_mad = c(F, T), resamp_when = "after_vs", which = "None"), resamp_methods)
resamp_methods$run <- 1:nrow(resamp_methods)
resamp_methods <- filter(resamp_methods, resamp_when == "before_vs", algo == "rf", which != "None")
resamp_methods <- resamp_methods[13,]

(cv_res <- ddply(resamp_methods, .(run), function(method) {
  cat("\n")
  print(method)

  loocv_res <- ddply(med2, .(id), function(r) {
    med_train <- subset(med2, id != r$id)

    if (method$do_mad) {
      train_mads <- apply(med_train[, med_preds], 2, function(g) mean(abs(g - median(g))))
      train_top_genes <- names(tail(sort(train_mads), 15000))
    } else {
      train_top_genes <- top_genes
    }

    if (method$resamp_when == "before_vs") {
      train_genes <- train_top_genes
    } else {
      train_genes <- top_genes_loo[[r$id]]
    }

    if (method$which != "None") {
      resamp_form <- as.formula(paste0("type ~ ", paste0(train_genes, collapse = "+")))
      response_vec <- as.factor(ifelse(med_train$type == "MLBCL", 1, 0))

      if (method$which == "Over2") {
        med_train_ub <- ubOver(med_train[, train_genes], response_vec, k = 2)
      } else if (method$which == "Over4") {
        med_train_ub <- ubOver(med_train[, train_genes], response_vec, k = 4)
      } else if (method$which == "Over0") {
        med_train_ub <- ubOver(med_train[, train_genes], response_vec, k = 0)
      } else if (method$which == "Under25") {
        med_train_ub <- ubUnder(med_train[, train_genes], response_vec, method = "percPos", perc = 25)
      } else if (method$which == "Under50") {
        med_train_ub <- ubUnder(med_train[, train_genes], response_vec, method = "percPos", perc = 50)
      } else if (method$which == "SMOTE_2_2") {
        med_train_ub <- ubSMOTE(med_train[, train_genes], response_vec, perc.over = 200, perc.under = 200)
      } else if (method$which == "SMOTE_3_1") {
        med_train_ub <- ubSMOTE(med_train[, train_genes], response_vec, perc.over = 300, perc.under = 100)
      } else if (method$which == "SMOTE_1_3") {
        med_train_ub <- ubSMOTE(med_train[, train_genes], response_vec, perc.over = 100, perc.under = 300)
      } else if (method$which == "Tomek") {
        med_train_ub <- ubTomek(med_train[, train_genes], response_vec)
      } else if (method$which == "NCL3") {
        med_train_ub <- ubNCL(med_train[, train_genes], response_vec, k = 3)
      } else if (method$which == "NCL4") {
        med_train_ub <- ubNCL(med_train[, train_genes], response_vec, k = 4)
      } else if (method$which == "SMOTE_2_2_Tomek") {
        med_train_ub <- ubSMOTE(med_train[, train_genes], response_vec, perc.over = 200, perc.under = 200)
        med_train_ub <- ubTomek(med_train_ub$X, med_train_ub$Y)
      } else if (method$which == "Tomek_SMOTE_2_2") {
        med_train_ub <- ubTomek(med_train[, train_genes], response_vec)
        med_train_ub <- ubSMOTE(med_train_ub$X, med_train_ub$Y, perc.over = 200, perc.under = 200)
      } else if (method$which == "SMOTE_2_2_NCL3") {
        med_train_ub <- ubNCL(med_train[, train_genes], response_vec, k = 3)
        med_train_ub <- ubSMOTE(med_train_ub$X, med_train_ub$Y, perc.over = 200, perc.under = 200)
      } else if (method$which == "NCL3_SMOTE_2_2") {
        med_train_ub <- ubSMOTE(med_train[, train_genes], response_vec, perc.over = 200, perc.under = 200)
        med_train_ub <- ubNCL(med_train_ub$X, med_train_ub$Y, k = 3)
      } else if (method$which == "ROSE") {
        med_train <- ROSE(resamp_form, med_train)$data
      } else if (method$which == "ADASYN") {
        med_train <- adasyn(resamp_form, med_train, beta = 0.5, do_parallel = F)
      } else if (method$which == "bSMOTE") {
        med_train <- bSMOTE(resamp_form, med_train, k = 5, s = 5, do_parallel = F)
      }

      if (!(method$which %in% c("ROSE", "ADASYN", "bSMOTE"))) {
        med_train <- cbind(type = med_train_ub$Y, med_train_ub$X)
        med_train$type <- as.factor(ifelse(med_train$type == "1", "MLBCL", "DLBCL"))
      }

      if (method$resamp_when == "before_vs") train_genes <- get_snr_genes(med_train, train_genes, m)
    }

    if (method$algo == "nb") {
      fold_res <- NaiveBayes(x = med_train[, train_genes], grouping = med_train$type, usekernel = T)
      r_nb_predict <- predict(fold_res, r[, train_genes])$class
    } else if (method$algo == "rf") {
      fold_res <- randomForest(x = as.matrix(med_train[, train_genes]), y = med_train$type, ntree = 501)
      r_nb_predict <- predict(fold_res, r[, train_genes])
    }

    return(data.frame(obs = r$type, pred = r_nb_predict))
  }, .parallel = T)

  (cm <- confusionMatrix(loocv_res$pred, loocv_res$obs, positive = "MLBCL"))

  f1 <- unname(2 * (cm$byClass["Pos Pred Value"] * cm$byClass["Sensitivity"]) / (cm$byClass["Pos Pred Value"] + cm$byClass["Sensitivity"]))
  g <- sqrt(cm$byClass["Sensitivity"] * cm$byClass["Specificity"])

  data.frame(algo = method$algo, do_mad = method$do_mad, resamp_when = method$resamp_when, which = method$which, tp = cm$table["MLBCL", "MLBCL"], fn = cm$table["MLBCL", "DLBCL"], fp = cm$table["DLBCL", "MLBCL"], tn = cm$table["DLBCL", "DLBCL"], accuracy = cm$overall["Accuracy"], bal_accuracy = cm$byClass["Balanced Accuracy"], f1 = f1, g = g, kappa = cm$overall["Kappa"], sensitivity = cm$byClass["Sensitivity"], specificity = cm$byClass["Specificity"], ppv = cm$byClass["Pos Pred Value"], npv = cm$byClass["Neg Pred Value"])
}, .inform = T, .progress = "text"))

save(cv_res_coll, file = "./cv_res.rdata")

cv_res_print <- cbind(cv_res[, 1:5], round(cv_res[, -5:-1], 2))[, c(-1, -3, -6:-9)]
print(xtable(cv_res_print), include.rownames = F)
