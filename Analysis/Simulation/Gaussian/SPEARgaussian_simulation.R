### Load Libraries:
# lasso
library(glmnet)
# MOFA
require(parallel)
library(reticulate)
library(MOFA2)
predict <- stats::predict

###########################################
#     REPRODUCIBILITY NOTE FOR USER:      #
#     CHANGE THE PYTHON PATH HERE:        #
python_path <- "/opt/anaconda3/bin/python"#
###########################################
use_python(python_path, required = TRUE)

# Path to save results to:
save.path <- "~/SPEAR_gaussian_results"

# SPEAR
library(SPEAR)

# initial seed:
seed <- 123
num_reps <- 1

# Functions:
preparation <- function(Y,  X, family, pattern_samples = NULL, pattern_assays = NULL,
                        path.type = "assay", bx = NULL,
                        other.path = NULL){
  ###check input types
  if(is.null(dim(Y))){
    Y = matrix(Y, ncol = 1)
  }
  py = ncol(Y)
  nclasses = rep(2, ncol(Y))
  for(j in 1:py){
    if(family != 0){
      labels = sort(unique(Y[,j]))
      labels.correct = 0:(length(labels)-1)
      if(sum(labels!= labels.correct) != 0){
        stop("class labels are not consecutive integers starting from 0")
      }
      nclasses[j] = length(labels)
    }else{
    }
  }
  if(!(path.type %in% c("assay", "none", "other"))){
    stop("Unsupported grouping structure.")
  }else if(path.type == "other" & is.null(other.path)){
    stop("other.grouping is missing.")
  }
  ##prepare the data
  X_ = X[[1]]
  n = nrow(X[[1]])
  p = rep(0, length(X))
  for(d in 1:length(X)){
    p[d] = ncol(X[[d]])
    if(d > 1){
      X_ = cbind(X_, X[[d]])
    }
  }
  bx_aggregated = NULL
  if(!is.null(bx)){
    bx_aggregated = bx[[1]]
    for(d in 1:length(X)){
      if(d > 1){
        bx_aggregated = cbind(bx_aggregated, bx[[d]])
      }
    }
    bx_aggregated = t(bx_aggregated)
  }
  functional_path = list()
  if(path.type == "assay"){
    for(i in 1:length(p)){
      if(i == 1){
        functional_path[[i]]  = 1:p[i]
      }else{
        functional_path[[i]] = (sum(p[1:(i-1)])+1):sum(p[1:(i)])
      }
      functional_path[[i]]  = functional_path[[i]]
    }
  }else if(group.type == "other"){
    functional_path= other.grouping;
  }else{
    functional_path[[1]] = c(1:sum(p));
  }
  px = ncol(X_)
  if(is.null(pattern_samples) | is.null(pattern_assays)){
    pattern_samples = list()
    pattern_assays = list()
    pattern_samples[[1]] = c(1:n)
    pattern_assays[[1]] = c(1:length(p))
  }else if(length(pattern_samples)!=length(pattern_assays)){
    stop("feature patterns and sample patterns do not match!")
  }else{
    tmp1 = pattern_samples[[1]]
    tmp2 = pattern_samples[[1]]
    for(k in 1:length(pattern_samples)){
      tmp1 = intersect(tmp1, pattern_samples[[k]])
      tmp2 = sort(union(tmp2, pattern_samples[[k]]))
      if(length(tmp1) > 0 | length(tmp2)!=n){
        stop("pattern_samples is not a partition of all samples!")
      }
    }
  }
  pattern_features = list()
  psum = cumsum(p)
  psum = c(0, psum)
  for(k in 1:length(pattern_assays)){
    pattern_features[[k]] = c(NA)
    for(l in pattern_assays[[k]]){
      pattern_features[[k]] = c(pattern_features[[k]], (psum[l]+1):psum[l+1])
    }
    pattern_features[[k]] = pattern_features[[k]][-1]
  }
  return(list(Y = Y, X = X_, functional_path = functional_path,
              pattern_samples = pattern_samples, pattern_features = pattern_features,
              nclasses = nclasses, bx_aggregated = bx_aggregated))
}
dataGen <- function(N = 500, Ntest = 2000, P = 500, D = 4, iteration = NULL, seed = 123, num_factors = 5, c = 1, 
                    pi = 0.2, eta = 1, num_specific =D-2, Ymodel = "factor", family = 0, factors.influencing.y = 2,
                    pi_reg = 0.05){
  
  set.seed(seed)
  Theta0 = list(); Gamma0 = list()
  X = list(); Xte = list()
  for(d in 1:D){
    Theta0[[d]] = matrix(rnorm(P*num_factors, sd = 1), ncol = P) * c
    Gamma0[[d]] = matrix(rbinom(P*num_factors, size = 1, prob  = pi), ncol = P)
  }
  if(num_specific != D){
    print("Assigning factors to datasets:")
    for(k in 1:num_factors){
      ii = sample(1:D, num_specific)
      print(paste0("Factor", k, ": ", paste(ii, collapse = ", ")))
      for(d in 1:D){
        if(!(d%in%ii)){
          Gamma0[[d]][k,] = 0
          Theta0[[d]] = Theta0[[d]] * Gamma0[[d]]
        }
      }
    }
  }
  U0 = matrix(rnorm(N*num_factors), nrow = N)
  U0te = matrix(rnorm(Ntest*num_factors), nrow = Ntest)
  X = list(); Xte = list()
  scale_mean.X = c(); scale_sd.X = c()
  for(d in 1:D){
    X[[d]] = scale(U0 %*% Theta0[[d]]+matrix(rnorm(N*P), ncol = P))
    Xte[[d]] = scale(U0te %*% Theta0[[d]]+matrix(rnorm(Ntest*P), ncol = P))
    tmp1 = apply(X[[d]], 2, mean);
    tmp2 =  apply(X[[d]], 2, sd);
    scale_mean.X =c(scale_mean.X, tmp1)
    scale_sd.X = c(scale_sd.X, tmp2)
    X[[d]] = t(apply(X[[d]], 1, function(z) (z - tmp1)/tmp2))
    Xte[[d]] = t(apply(Xte[[d]], 1, function(z) (z - tmp1)/tmp2))
  }
  if(Ymodel == "factor"){
    Y = rowSums(U0[,1:factors.influencing.y, drop=FALSE])*sqrt(eta/factors.influencing.y) + rnorm(N)
    Yte = rowSums(U0te[,1:factors.influencing.y, drop=FALSE])* sqrt(eta/factors.influencing.y)+rnorm(Ntest)
    Ytruth = rowSums(U0[,1:factors.influencing.y, drop=FALSE])*sqrt(eta/factors.influencing.y)
    Ytruth_te =  rowSums(U0te[,1:factors.influencing.y, drop=FALSE])* sqrt(eta/factors.influencing.y)
  }else{
    Xcombine = X[[1]]
    Xcombine.te = Xte[[1]]
    for(d in 2:D){
      Xcombine = cbind(Xcombine, X[[d]])
      Xcombine.te = cbind(Xcombine.te, Xte[[d]])
    }
    beta =rnorm(ncol(Xcombine)) * rbinom(ncol(Xcombine), size = 1, prob = pi_reg) 
    beta = beta/sqrt(sum(beta^2)) * sqrt(eta)
    Ytruth =  Xcombine%*%beta
    Ytruth_te =  Xcombine.te%*% beta
    Y = Ytruth + rnorm(N)
    Yte = Ytruth_te+rnorm(Ntest)
  }
  scale_mean.y = mean(Y); scale_sd.y = sd(Y);
  y = (Y - scale_mean.y)/scale_sd.y; yte = (Yte - scale_mean.y)/scale_sd.y;
  ytruth = (Ytruth - scale_mean.y)/scale_sd.y; ytruth.te = (Ytruth_te - scale_mean.y)/scale_sd.y;
  data.tr = preparation(Y = y, X = X, family = family, path.type = "assay")
  data.te = preparation(Y = yte, X = Xte, family = family, path.type = "assay")
  data.tr$truth = ytruth
  data.tr$xlist = X
  data.tr$U = U0
  data.te$truth = ytruth.te
  data.te$U = U0te
  data.te$xlist = Xte
  return(list(data.tr = data.tr, data.te = data.te))
}

# Get arguments:
args = commandArgs(trailingOnly=TRUE)
print(length(args))
if(length(args) == 0){
  # Data simulation parameters:
  K = 5 # number of factors
  D = 4 # number of omics datasets to simulate
  N = 100
  P = 100
  c = .3
  c.original = c
  c = sqrt(c *log(P*D)/N) 
  Ymodel = "specificfactors" # parameter indicating that Y (gaussian response) is generated from the factors
  Ymodel.original = Ymodel
  Ntest = 2000 # number of subjects in the testing set
  factors.influencing.y = 2 # how many of the K factors should be used to influence the response?
  family = 0 # When response is gaussian, use family = 0. For ordinal response, use family = 1
  pi = .2 # amount of sparsity for features (P(feature having signal) = pi)
  
  
  if(Ymodel == "XY"){
    factors.influencing.y = 2
    num_specific = 1
  } else if(Ymodel == "sharedfactors"){
    num_specific = D
    Ymodel = "factor"
  } else {
    num_specific = 2
    Ymodel = "factor"
  }
} else {
  # Data simulation parameters:
  K = 5 # number of factors
  D = 4 # number of omics datasets to simulate
  N = as.numeric(args[1]) # number of subjects in the training set
  P = as.numeric(args[2]) # number of features to simulate per omics dataset
  c = as.numeric(args[3]) # signal-to-noise ratio. 1 is low, 5 is high
  c.original = c
  c = sqrt(c *log(P*D)/N) 
  Ymodel = as.character(args[4]) # parameter indicating that Y (gaussian response) is generated from the factors
  Ymodel.original = Ymodel
  Ntest = 2000 # number of subjects in the testing set
  factors.influencing.y = 2 # how many of the K factors should be used to influence the response?
  family = 0 # When response is gaussian, use family = 0. For ordinal response, use family = 1
  pi = .2 # amount of sparsity for features (P(feature having signal) = pi)
  
  
  if(Ymodel == "XY"){
    factors.influencing.y = 2
    num_specific = 1
  } else if(Ymodel == "sharedfactors"){
    num_specific = D
    Ymodel = "factor"
  } else {
    num_specific = 2
    Ymodel = "factor"
  }
}

final.res <- list()
for(iteration in 1:num_reps){
  cat(paste0("+++++++++++++++++++++++\nStarting Iteration", iteration, ":\n+++++++++++++++++++++++\n"))
  seed <- seed + iteration
  
  # Data Generation:
  sim.data = dataGen(N = N, 
                     Ntest = Ntest, 
                     P = P, 
                     D = D, 
                     iteration = NULL, 
                     seed = seed, 
                     num_factors = K, 
                     c = c, 
                     pi = pi, 
                     eta = 1, 
                     num_specific =num_specific, 
                     Ymodel = Ymodel,
                     family = 0,
                     factors.influencing.y = factors.influencing.y,
                     pi_reg = N/(P*D*log(P * D) * log(P * D)))
  
  family <- "gaussian"
  X <- sim.data$data.tr$xlist
  Y <- sim.data$data.tr$Y
  X.te <- sim.data$data.te$xlist
  Y.te <- sim.data$data.te$Y
  set.seed(seed + 100)
  fold.ids <- sample(1:5, nrow(Y), replace = TRUE)
  set.seed(seed)
  
  total.results <- list()
  
  # Lasso:
  lasso_fit = cv.glmnet(x = do.call("cbind", X), y = Y, foldid = fold.ids)
  lasso.preds.tr = predict(lasso_fit, do.call("cbind", X), s = "lambda.min")
  lasso.preds.te = predict(lasso_fit, do.call("cbind", X.te), s = "lambda.min")
  total.results$lasso <- list(train = lasso.preds.tr,
                              test = lasso.preds.te)
  
  
  
  
  # SPEAR:
  num.factors.spear <- 5
  SPEAR_fit = SPEAR::make.SPEARobject(X = X,
                                      Y = Y,
                                      family = family,
                                      num_factors = num.factors.spear)
  
  
  SPEAR_fit$add.data(X = X.te, Y = Y.te, name = "test")
  SPEAR_fit$run.cv.spear(fold.ids = fold.ids)$cv.evaluate()
  SPEAR.preds <- list()
  for(w in 1:nrow(SPEAR_fit$params$weights)){
    SPEAR_fit$set.weights(w.idx = w)
    SPEAR.preds[[paste0("widx", w)]] <- list(train = SPEAR_fit$get.predictions(cv = FALSE),
                                             cv = SPEAR_fit$get.predictions(cv = TRUE),
                                             test = SPEAR_fit$get.predictions(data = "test"))
  }
  total.results$SPEAR <- SPEAR.preds
  
  
  
  # MOFA:
  x.mofa <- list()
  for(d in 1:length(X)){
    x.mofa[[d]] <- t(X[[d]])
  }
  MOFAobject <- MOFA2::create_mofa(x.mofa)
  # MOFA+ specific parameters:
  data_opts <- get_default_data_options(MOFAobject)
  model_opts <- get_default_model_options(MOFAobject)
  model_opts$num_factors <- 10
  train_opts <- get_default_training_options(MOFAobject)
  train_opts$seed <- seed
  MOFAobject <- prepare_mofa(
    object = MOFAobject,
    data_options = data_opts,
    model_options = model_opts,
    training_options = train_opts
  )
  # Running MOFA+
  MOFA_res <- run_mofa(MOFAobject, save_data = FALSE)
  if(ncol(MOFA_res@expectations$Z[[1]]) < 2){
    MOFA.factors <- cbind(MOFA_res@expectations$Z[[1]], MOFA_res@expectations$Z[[1]])
    colnames(MOFA.factors) <- c("Factor1", "Factor2")
  } else {
    MOFA.factors <- MOFA_res@expectations$Z[[1]]
  }
  MOFA.coefs = matrix(NA, ncol = ncol(MOFA.factors), nrow = ncol(do.call("cbind", X)))
  for(k in 1:ncol(MOFA.factors)){
    tmp =  cv.glmnet(x = do.call("cbind", X), y = MOFA.factors[,k])
    MOFA.coefs[,k] = coef(tmp, s = "lambda.min")[-1]
  }
  MOFA.factors.tr = do.call("cbind", X) %*% MOFA.coefs
  MOFA.factors.te = do.call("cbind", X.te) %*% MOFA.coefs
  MOFA_lasso_res = cv.glmnet(MOFA.factors.tr, Y)
  MOFA.train = predict(MOFA_lasso_res, MOFA.factors.tr, s = "lambda.min")
  MOFA.test = predict(MOFA_lasso_res, MOFA.factors.te, s = "lambda.min")
  total.results$MOFA <- list(train = MOFA.train, cv = NA, test = MOFA.test)
  
  
  
  ### Calculate results:
  error.res <- list()
  for(i in 1:length(total.results)){
    model <- total.results[[i]]
    if(names(total.results)[i] == "SPEAR"){
      temp <- list()
      for(w in 1:length(model)){
        train.err <- as.numeric(mean((total.results$SPEAR[[w]]$train - Y)^2))
        train.cv.err <- as.numeric(mean((total.results$SPEAR[[w]]$cv - Y)^2))
        test.err <- as.numeric(mean((total.results$SPEAR[[w]]$test - Y.te)^2))
        temp[[paste0(names(total.results)[i], "_widx", w)]] <- list(train = train.err, cv = train.cv.err, test = test.err, w.idx = w)
      }
      error.df <- do.call("rbind", temp)
    } else {
      train.err <- mean((model$train - Y)^2)
      test.err <- mean((model$test - Y.te)^2)
      error.res[[names(total.results)[i]]] <- list(train = train.err, cv = NA, test = test.err, w.idx = NA)
    }
  }
  
  # Correlations:
  cat("\nDOING CORRELATIONS... \n")
  MOFA.factors <- MOFA2::get_factors(MOFA_res)$group1
  colnames(MOFA.factors) <- paste0("MOFA.", colnames(MOFA.factors))
  total.factors <- MOFA.factors
  if(Ymodel != "XY"){
    U.true <- sim.data$data.tr$U
    colnames(U.true) <- paste0("True.Factor", 1:ncol(U.true))
    rownames(U.true) <- rownames(total.factors)
    total.factors <- cbind(total.factors, U.true)
  }
  for(i in 1:5){
    SPEAR.factors <- SPEAR_fit$set.weights(w.idx = i)$get.factor.scores()
    colnames(SPEAR.factors) <- paste0("SPEAR.widx", i, ".", colnames(SPEAR.factors))
    rownames(SPEAR.factors) <- rownames(total.factors)
    total.factors <- cbind(total.factors, SPEAR.factors)
  }
  total.cor.mat <- matrix(0, nrow = ncol(total.factors), ncol = ncol(total.factors))
  colnames(total.cor.mat) <- rownames(total.cor.mat) <- colnames(total.factors)
  total.cor.mat.pvals <- matrix(0, nrow = ncol(total.factors), ncol = ncol(total.factors))
  colnames(total.cor.mat.pvals) <- rownames(total.cor.mat.pvals) <- colnames(total.factors)
  for(i in 1:nrow(total.cor.mat)){
    for(j in 1:ncol(total.cor.mat)){
      x = colnames(total.factors)[i]
      y = colnames(total.factors)[j]
      cor.res <- cor.test(total.factors[,x], total.factors[,y], method = "spearman")
      total.cor.mat[i,j] <- as.numeric(cor.res$estimate)
      total.cor.mat.pvals[i,j] <- as.numeric(cor.res$p.value)
    }
  }
  
  
  # Simplify results:
  for(i in 1:length(error.res)){
    error.df <- rbind(error.df, c(error.res[[i]]$train,
                                  error.res[[i]]$cv,
                                  error.res[[i]]$test,
                                  error.res[[i]]$w.idx))
    rownames(error.df)[nrow(error.df)] <- names(error.res)[i]
  }
  error.df <- as.data.frame(error.df)
  
  final.res[[paste0("Iteration", iteration)]] <- list(errors = error.df, correlations = list(cors = total.cor.mat, p.values = total.cor.mat.pvals))
  
}

saveRDS(final.res, paste0(save.path, "/c", c.original, "_", Ymodel.original, ".rds"))


