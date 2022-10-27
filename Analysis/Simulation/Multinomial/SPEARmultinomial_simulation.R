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
save.path <- "~/SPEAR_multinomial_results"

# SPEAR
library(SPEAR)

# initial seed:
seed <- 123
num_reps <- 10
K <- 5
D <- 4
N.test = 2000
num.factors.spear <- 5
eta <- 10
family = "multinomial"

multi.centers.x = c(-3, -3, -3, -1, -1, -1, 1, 1, 1)/2
multi.centers.y = c(-3, -1, 1, -3, -1, 1, -3, -1, 1)/2


simulate.data <- function(N = 100, D = 3, P = 100, c1 = 3, c2 = 1, sparsity = .2,
                          method = "factor", num.factors = 5, family = "gaussian", factors.influencing.X = 2, factors.influencing.Y = 2, ordinal.centers = c(-1.2,-1, 0, 1,1.2), multi.centers.x = c(-1, -1, 1, 1)/sqrt(2), multi.centers.y = c(-1, 1, -1, 1)/sqrt(2),
                          N.test = 1000, seed = 123){
  # Set the seed:
  set.seed(seed)
  
  if(family == "gaussian"){
    Theta0 = list(); Gamma0 = list()
    X = list(); Xte = list()
    for(d in 1:D){
      Theta0[[d]] = matrix(rnorm(P*num.factors, sd = 1), ncol = P) * c2
      Gamma0[[d]] = matrix(rbinom(P*num.factors, size = 1, prob  = sparsity), ncol = P)
    }
    if(factors.influencing.X != D){
      for(k in 1:num.factors){
        ii = sample(1:D, factors.influencing.X)
        for(d in 1:D){
          if(!(d%in%ii)){
            Gamma0[[d]][k,] = 0
            Theta0[[d]] = Theta0[[d]] * Gamma0[[d]]
          }
        }
      }
    }
    U0 = matrix(rnorm(N*num.factors), nrow = N)
    U0te = matrix(rnorm(N.test*num.factors), nrow = N.test)
    X = list(); Xte = list()
    scale_mean.X = c(); scale_sd.X = c()
    for(d in 1:D){
      X[[d]] = scale(U0 %*% Theta0[[d]]+matrix(rnorm(N*P), ncol = P))
      Xte[[d]] = scale(U0te %*% Theta0[[d]]+matrix(rnorm(N.test*P), ncol = P))
      tmp1 = apply(X[[d]], 2, mean);
      tmp2 =  apply(X[[d]], 2, sd);
      scale_mean.X =c(scale_mean.X, tmp1)
      scale_sd.X = c(scale_sd.X, tmp2)
      X[[d]] = t(apply(X[[d]], 1, function(z) (z - tmp1)/tmp2))
      Xte[[d]] = t(apply(Xte[[d]], 1, function(z) (z - tmp1)/tmp2))
    }
    if(method == "factor"){
      Y = rowSums(U0[,1:factors.influencing.Y, drop=FALSE])*sqrt(c2/factors.influencing.Y) + rnorm(N)
      Yte = rowSums(U0te[,1:factors.influencing.Y, drop=FALSE])* sqrt(c2/factors.influencing.Y)+rnorm(N.test)
      Ytruth = rowSums(U0[,1:factors.influencing.Y, drop=FALSE])*sqrt(c2/factors.influencing.Y)
      Ytruth_te =  rowSums(U0te[,1:factors.influencing.Y, drop=FALSE])* sqrt(c2/factors.influencing.Y)
    }else{
      Xcombine = X[[1]]
      Xcombine.te = Xte[[1]]
      for(d in 2:D){
        Xcombine = cbind(Xcombine, X[[d]])
        Xcombine.te = cbind(Xcombine.te, Xte[[d]])
      }
      bc2 =rnorm(ncol(Xcombine)) * rbinom(ncol(Xcombine), size = 1, prob = .05) 
      bc2 = bc2/sqrt(sum(bc2^2)) * sqrt(c2)
      Ytruth =  Xcombine%*%bc2
      Ytruth_te =  Xcombine.te%*% bc2
      Y = Ytruth + rnorm(N)
      Yte = Ytruth_te+rnorm(N.test)
    }
    scale_mean.y = mean(Y); scale_sd.y = sd(Y);
    y = (Y - scale_mean.y)/scale_sd.y; yte = (Yte - scale_mean.y)/scale_sd.y;
    ytruth = (Ytruth - scale_mean.y)/scale_sd.y; ytruth.te = (Ytruth_te - scale_mean.y)/scale_sd.y;
    data.tr = preparation(Y = y, X = X, family = 0, path.type = "assay")
    data.te = preparation(Y = yte, X = Xte, family = 0, path.type = "assay")
    data.tr$truth = ytruth
    data.tr$Xlist = X
    data.tr$U = U0
    data.te$truth = ytruth.te
    data.te$U = U0te
    data.te$Xlist = Xte
    
    return(list(train = data.tr, test = data.te))
    
  } else if(family == "ordinal"){
    Theta0 = list(); Gamma0 = list()
    X = list(); Xte = list()
    for(d in 1:D){
      Theta0[[d]] = matrix(rnorm(P*num.factors, sd = 1), ncol = P) * c1
      Gamma0[[d]] = matrix(rbinom(P*num.factors, size = 1, prob  = sparsity), ncol = P)
    }
    for(k in 1:num.factors){
      ii = sample(1:D, factors.influencing.X)
      for(d in 1:D){
        if(!(d%in%ii)){
          Gamma0[[d]][k,] = 0
          Theta0[[d]] = Theta0[[d]] * Gamma0[[d]]
        }
      }
    }
    U0 = matrix(rnorm(N*num.factors), nrow = N)
    U0te = matrix(rnorm(N.test*num.factors), nrow = N.test)
    ###create Y based on specified signals
    nclasses <- length(ordinal.centers)
    Y = sample(0:(nclasses-1), N, replace = T)
    Yte = sample(0:(nclasses-1), N.test, replace = T)
    mu = ordinal.centers[Y+1]
    mu.te = ordinal.centers[Yte+1]
    ##add noise
    for(i in 1:factors.influencing.Y){
      U0[,i] = mu+rnorm(N) * sqrt(1/(c2+1))
      U0te[,i] = mu.te+rnorm(N.test) * sqrt(1/(c2+1))
    }
    ms = apply(U0,2,mean)
    sds = apply(U0,2,sd)
    U0 = t(apply(U0, 1, function(z) (z-ms)/sds))
    U0te = t(apply(U0te, 1, function(z) (z-ms)/sds))
    mu = (mu-ms[1])/sds[1]
    mu.te = (mu.te-ms[1])/sds[1]
    X = list(); Xte = list()
    scale_mean.X = c(); scale_sd.X = c()
    for(d in 1:D){
      X[[d]] = scale(U0 %*% Theta0[[d]]+matrix(rnorm(N*P), ncol = P))
      Xte[[d]] = scale(U0te %*% Theta0[[d]]+matrix(rnorm(N.test*P), ncol = P))
      tmp1 = apply(X[[d]], 2, mean);
      tmp2 =  apply(X[[d]], 2, sd);
      scale_mean.X =c(scale_mean.X, tmp1)
      scale_sd.X = c(scale_sd.X, tmp2)
      X[[d]] = t(apply(X[[d]], 1, function(z) (z - tmp1)/tmp2))
      Xte[[d]] = t(apply(Xte[[d]], 1, function(z) (z - tmp1)/tmp2))
    }
    data.tr = preparation(Y = Y, X = X, family = 2, path.type = "assay")
    data.te = preparation(Y = Yte, X = Xte, family = 2,path.type = "assay")
    data.tr$truth = mu
    data.tr$Xlist = X
    data.tr$U = U0
    data.te$truth = mu.te
    data.te$U = U0te
    data.te$Xlist = Xte
    data.tr$centers = ordinal.centers
    data.te$centers = ordinal.centers
    Y.multi <- reshape2::dcast(data.frame(x = 1:nrow(data.tr$Y), y = data.tr$Y), x ~ y, length, value.var = "y")
    data.tr$Y.onehot <- Y.multi[,2:ncol(Y.multi)]
    Yte.multi <- reshape2::dcast(data.frame(x = 1:nrow(data.te$Y), y = data.te$Y), x ~ y, length, value.var = "y")
    data.te$Y.onehot <- Yte.multi[,2:ncol(Yte.multi)]
    
    return(list(train = data.tr, test = data.te))
    
  } else if(family == "multinomial"){
    Theta0 = list(); Gamma0 = list()
    X = list(); Xte = list()
    for(d in 1:D){
      Theta0[[d]] = matrix(rnorm(P*num.factors, sd = 1), ncol = P) * c1
      Gamma0[[d]] = matrix(rbinom(P*num.factors, size = 1, prob  = sparsity), ncol = P)
    }
    for(k in 1:num.factors){
      ii = sample(1:D, factors.influencing.X)
      for(d in 1:D){
        if(!(d%in%ii)){
          Gamma0[[d]][k,] = 0
          Theta0[[d]] = Theta0[[d]] * Gamma0[[d]]
        }
      }
    }
    U0 = matrix(rnorm(N*num.factors), nrow = N)
    U0te = matrix(rnorm(N.test*num.factors), nrow = N.test)
    ###create Y based on specified signals
    nclasses <- length(multi.centers.x)
    Y = sample(0:(nclasses-1), N, replace = T)
    Yte = sample(0:(nclasses-1), N.test, replace = T)
    mu = matrix(0, N, 2)
    mu.te = matrix(0, N.test, 2)
    mu[,1] = multi.centers.x[Y+1];  mu[,2] = multi.centers.y[Y+1]
    mu.te[,1] = multi.centers.x[Yte+1];  mu.te[,2] = multi.centers.y[Yte+1]
    ##add noise
    U0[,c(1:2)] = mu+matrix(rnorm(N*2),ncol=2) * sqrt(1/(c2+1))
    U0te[,c(1:2)] = mu.te+matrix(rnorm(N.test*2),ncol=2) * sqrt(1/(c2+1))
    ms = apply(U0,2,mean)
    sds = apply(U0,2,sd)
    U0 = t(apply(U0, 1, function(z) (z-ms)/sds))
    U0te = t(apply(U0te, 1, function(z) (z-ms)/sds))
    mu = (mu-ms[1])/sds[1]
    mu.te = (mu.te-ms[1])/sds[1]
    X = list(); Xte = list()
    scale_mean.X = c(); scale_sd.X = c()
    for(d in 1:D){
      X[[d]] = scale(U0 %*% Theta0[[d]]+matrix(rnorm(N*P), ncol = P))
      Xte[[d]] = scale(U0te %*% Theta0[[d]]+matrix(rnorm(N.test*P), ncol = P))
      tmp1 = apply(X[[d]], 2, mean);
      tmp2 =  apply(X[[d]], 2, sd);
      scale_mean.X =c(scale_mean.X, tmp1)
      scale_sd.X = c(scale_sd.X, tmp2)
      X[[d]] = t(apply(X[[d]], 1, function(z) (z - tmp1)/tmp2))
      Xte[[d]] = t(apply(Xte[[d]], 1, function(z) (z - tmp1)/tmp2))
    }
    data.tr = preparation(Y = Y, X = X, family = 2, path.type = "assay")
    data.te = preparation(Y = Yte, X = Xte, family = 2,path.type = "assay")
    data.tr$truth = mu
    data.tr$Xlist = X
    data.tr$U = U0
    data.te$truth = mu.te
    data.te$U = U0te
    data.te$Xlist = Xte
    data.tr$centers.x = multi.centers.x
    data.tr$centers.y = multi.centers.y
    data.te$centers.x = multi.centers.x
    data.te$centers.y = multi.centers.y
    Y.multi <- reshape2::dcast(data.frame(x = 1:nrow(data.tr$Y), y = data.tr$Y), x ~ y, length, value.var = "y")
    data.tr$Y.onehot <- Y.multi[,2:ncol(Y.multi)]
    Yte.multi <- reshape2::dcast(data.frame(x = 1:nrow(data.te$Y), y = data.te$Y), x ~ y, length, value.var = "y")
    data.te$Y.onehot <- Yte.multi[,2:ncol(Yte.multi)]
    
    return(list(train = data.tr, test = data.te))
    
  } else {
    stop("ERROR: family parameter provided (", family, ") not recognized. Can be 'gaussian', 'ordinal', or 'multinomial'.")
  }
}


# Get arguments:
args = commandArgs(trailingOnly=TRUE)
print(length(args))
if(length(args) == 0){
  # Data simulation parameters:
  N = 500
  P = 500
  c = 1 # 3, 10, 20
  c.original = c
  c = sqrt(c *log(P*D)/N) 
  
} else {
  # Data simulation parameters:
  N = as.numeric(args[1]) # number of subjects in the training set
  P = as.numeric(args[2]) # number of features to simulate per omics dataset
  c = as.numeric(args[3]) # signal-to-noise ratio. 1 is low, 5 is high
  c.original = c
  c = sqrt(c *log(P*D)/N) 
}



final.res <- list()
for(iteration in 1:num_reps){
      cat(paste0("+++++++++++++++++++++++\nStarting Iteration", iteration, ":\n+++++++++++++++++++++++\n"))
      seed <- seed + iteration
      
      # Generate Data:
      sim.data <- simulate.data(family = family,
                                c2 = eta,
                                c1 = c,
                                multi.centers.x = multi.centers.x,
                                multi.centers.y = multi.centers.y,
                                N = N,
                                P = P,
                                D = D,
                                N.test = N.test,
                                seed = seed)
      
    plot(sim.data$train$U[,2] ~ sim.data$train$U[,1], col = sim.data$train$Y+1)
    plot(sim.data$test$U[,2] ~ sim.data$test$U[,1], col = sim.data$test$Y+1)
    
    # Make data accessible:    
    X <- sim.data$train$Xlist
    Y <- sim.data$train$Y
    Y.multi <- sim.data$train$Y.onehot
    X.te <- sim.data$test$Xlist
    Y.te <- sim.data$test$Y
    Y.te.multi <- sim.data$test$Y.onehot
    
    # Make SPEARobjects for gaussian and multinomial:
    SPEARgaussian <- make.SPEARobject(X = X,
                                      Y = Y,
                                      family = "gaussian",
                                      num_factors = 10,
                                      seed = seed)
    SPEARgaussian$add.data(X = X.te, Y = Y.te, name = "test")
    
    SPEARmulti <- make.SPEARobject(X = X,
                                   Y = Y.multi,
                                   family = "multinomial",
                                   num_factors = 10,
                                   seed = seed)
    SPEARmulti$add.data(X = X.te, Y = Y.te.multi, name = "test")
    
    
    # Make fold.ids:
    SPEARmulti$params$num.folds <- 5
    fold.ids <- SPEARmulti$generate.fold.ids(method = "balanced")
    
    ### RUN MODELS:

    # Run CV SPEARmulti and SPEARgaussian:
    SPEARgaussian$run.cv.spear(fold.ids = fold.ids)$cv.evaluate()
    SPEARmulti$run.cv.spear(fold.ids = fold.ids)$cv.evaluate()
    
    # Run Lasso gaussian and multinomial:
    lasso.gaus_fit = cv.glmnet(x = do.call("cbind", X), y = Y, foldid = fold.ids)
    lasso.multi_fit = cv.glmnet(x = do.call("cbind", X), y = as.matrix(Y.multi), foldid = fold.ids, family = "multinomial")

    # Run MOFA:
    # MOFA:
    train.data <- list()
    train.data$X <- sim.data$train$Xlist
    train.data$Y <- sim.data$train$Y.onehot
    test.data <- list()
    test.data$X <- sim.data$test$Xlist
    test.data$Y <- sim.data$test$Y.onehot
    
    x.mofa <- list()
    for(d in 1:length(train.data$X)){
      x.mofa[[d]] <- t(train.data$X[[d]])
    }
    MOFAobject <- MOFA2::create_mofa(x.mofa)
    # MOFA+ specific parameters:
    data_opts <- get_default_data_options(MOFAobject)
    model_opts <- get_default_model_options(MOFAobject)
    model_opts$num_factors <- 20
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
    MOFA.coefs = matrix(NA, ncol = ncol(MOFA.factors), nrow = ncol(do.call("cbind", train.data$X)))
    for(k in 1:ncol(MOFA.factors)){
      tmp =  glmnet::cv.glmnet(x = do.call("cbind", train.data$X), y = MOFA.factors[,k])
      MOFA.coefs[,k] = coef(tmp, s = "lambda.min")[-1]
    }
    MOFA.factors.tr = do.call("cbind", train.data$X) %*% MOFA.coefs
    MOFA.factors.te = do.call("cbind", test.data$X) %*% MOFA.coefs
    
    
    ### GET PREDICTIONS:    
    # SPEARgaussian
    # min
    SPEARgaussian$set.weights(method = "min")
    SPEAR.min.gaus.preds.tr <- SPEARgaussian$get.predictions(data = "train", cv = FALSE) + SPEARgaussian$fit$intercepts.scaled[[1]][SPEARgaussian$options$current.weight.idx,]
    SPEAR.min.gaus.preds.tr <- round(SPEAR.min.gaus.preds.tr, 0)
    SPEAR.min.gaus.preds.tr[SPEAR.min.gaus.preds.tr < 0] <- 0
    SPEAR.min.gaus.preds.tr[SPEAR.min.gaus.preds.tr > (length(multi.centers.x)-1)] <- length(multi.centers.x)-1
    SPEAR.min.gaus.preds.te <- SPEARgaussian$get.predictions(data = "test")  + SPEARgaussian$fit$intercepts.scaled[[1]][SPEARgaussian$options$current.weight.idx,]
    SPEAR.min.gaus.preds.te <- round(SPEAR.min.gaus.preds.te, 0)
    SPEAR.min.gaus.preds.te[SPEAR.min.gaus.preds.te < 0] <- 0
    SPEAR.min.gaus.preds.te[SPEAR.min.gaus.preds.te > (length(multi.centers.x)-1)] <- length(multi.centers.x)-1
    # sd
    SPEARgaussian$set.weights(method = "sd")
    SPEAR.sd.gaus.preds.tr <- SPEARgaussian$get.predictions(data = "train", cv = FALSE) + SPEARgaussian$fit$intercepts.scaled[[1]][SPEARgaussian$options$current.weight.idx,]
    SPEAR.sd.gaus.preds.tr <- round(SPEAR.sd.gaus.preds.tr, 0)
    SPEAR.sd.gaus.preds.tr[SPEAR.sd.gaus.preds.tr < 0] <- 0
    SPEAR.sd.gaus.preds.tr[SPEAR.sd.gaus.preds.tr > (length(multi.centers.x)-1)] <- length(multi.centers.x)-1
    SPEAR.sd.gaus.preds.te <- SPEARgaussian$get.predictions(data = "test")  + SPEARgaussian$fit$intercepts.scaled[[1]][SPEARgaussian$options$current.weight.idx,]
    SPEAR.sd.gaus.preds.te <- round(SPEAR.sd.gaus.preds.te, 0)
    SPEAR.sd.gaus.preds.te[SPEAR.sd.gaus.preds.te < 0] <- 0
    SPEAR.sd.gaus.preds.te[SPEAR.sd.gaus.preds.te > (length(multi.centers.x)-1)] <- length(multi.centers.x)-1
    # vanilla
    SPEARgaussian$set.weights(method = "vanilla")
    SPEAR.van.gaus.preds.tr <- SPEARgaussian$get.predictions(data = "train", cv = FALSE) + SPEARgaussian$fit$intercepts.scaled[[1]][SPEARgaussian$options$current.weight.idx,]
    SPEAR.van.gaus.preds.tr <- round(SPEAR.van.gaus.preds.tr, 0)
    SPEAR.van.gaus.preds.tr[SPEAR.van.gaus.preds.tr < 0] <- 0
    SPEAR.van.gaus.preds.tr[SPEAR.van.gaus.preds.tr > (length(multi.centers.x)-1)] <- length(multi.centers.x)-1
    SPEAR.van.gaus.preds.te <- SPEARgaussian$get.predictions(data = "test")  + SPEARgaussian$fit$intercepts.scaled[[1]][SPEARgaussian$options$current.weight.idx,]
    SPEAR.van.gaus.preds.te <- round(SPEAR.van.gaus.preds.te, 0)
    SPEAR.van.gaus.preds.te[SPEAR.van.gaus.preds.te < 0] <- 0
    SPEAR.van.gaus.preds.te[SPEAR.van.gaus.preds.te > (length(multi.centers.x)-1)] <- length(multi.centers.x)-1
    
    # SPEARmultinomial
    # min
    SPEARmulti$set.weights(method = "min")
    SPEAR.min.multi.preds.tr <- SPEARmulti$get.predictions(data = "train", cv = FALSE)$class.predictions
    SPEAR.min.multi.preds.te <- SPEARmulti$get.predictions(data = "test")$class.predictions
    # min
    SPEARmulti$set.weights(method = "sd")
    SPEAR.sd.multi.preds.tr <- SPEARmulti$get.predictions(data = "train", cv = FALSE)$class.predictions
    SPEAR.sd.multi.preds.te <- SPEARmulti$get.predictions(data = "test")$class.predictions
    # vanilla
    SPEARmulti$set.weights(method = "vanilla")
    SPEAR.van.multi.preds.tr <- SPEARmulti$get.predictions(data = "train", cv = FALSE)$class.predictions
    SPEAR.van.multi.preds.te <- SPEARmulti$get.predictions(data = "test")$class.predictions
    
    # Lasso gaussian
    lasso.gaus.preds.tr = predict(lasso.gaus_fit, do.call("cbind", X), s = "lambda.min")
    lasso.gaus.preds.tr <- round(lasso.gaus.preds.tr, 0)
    lasso.gaus.preds.tr[lasso.gaus.preds.tr < 0] <- 0
    lasso.gaus.preds.tr[lasso.gaus.preds.tr > (length(multi.centers.x)-1)] <- length(multi.centers.x)-1
    lasso.gaus.preds.te = predict(lasso.gaus_fit, do.call("cbind", X.te), s = "lambda.min")
    lasso.gaus.preds.te <- round(lasso.gaus.preds.te, 0)
    lasso.gaus.preds.te[lasso.gaus.preds.te < 0] <- 0
    lasso.gaus.preds.te[lasso.gaus.preds.te > (length(multi.centers.x)-1)] <- length(multi.centers.x)-1
    
    # Lasso multinomial
    lasso.multi.preds.tr = apply(predict(lasso.multi_fit, do.call("cbind", X), s = "lambda.min"), 1, which.max)-1
    lasso.multi.preds.te = apply(predict(lasso.multi_fit, do.call("cbind", X.te), s = "lambda.min"), 1, which.max)-1
    
    # MOFA:
    # Gaussian:
    MOFA.gaus_res = glmnet::cv.glmnet(MOFA.factors.tr, Y, foldid = fold.ids)
    MOFA.gaus.preds.tr = predict(MOFA.gaus_res, MOFA.factors.tr, s = "lambda.min")
    MOFA.gaus.preds.tr <- round(MOFA.gaus.preds.tr, 0)
    MOFA.gaus.preds.tr[MOFA.gaus.preds.tr < 0] <- 0
    MOFA.gaus.preds.tr[MOFA.gaus.preds.tr > (length(multi.centers.x)-1)] <- length(multi.centers.x)-1
    MOFA.gaus.preds.te = predict(MOFA.gaus_res, MOFA.factors.te, s = "lambda.min")
    MOFA.gaus.preds.te <- round(MOFA.gaus.preds.te, 0)
    MOFA.gaus.preds.te[MOFA.gaus.preds.te < 0] <- 0
    MOFA.gaus.preds.te[MOFA.gaus.preds.te > (length(multi.centers.x)-1)] <- length(multi.centers.x)-1
    # Multinomial:
    MOFA.multi_res = glmnet::cv.glmnet(MOFA.factors.tr, as.matrix(Y.multi), foldid = fold.ids, family = "multinomial")
    MOFA.multi.preds.tr = apply(predict(MOFA.multi_res, MOFA.factors.tr, s = "lambda.min"), 1, which.max)-1
    MOFA.multi.preds.te = apply(predict(MOFA.multi_res, MOFA.factors.te, s = "lambda.min"), 1, which.max)-1
    
    # Compile results:
    pred.df.tr <- data.frame(
      SPEAR.min.gaussian = SPEAR.min.gaus.preds.tr,
      SPEAR.sd.gaussian = SPEAR.sd.gaus.preds.tr,
      SPEAR.van.gaussian = SPEAR.van.gaus.preds.tr,
      MOFA.gaussian = MOFA.gaus.preds.tr,
      Lasso.gaussian = lasso.gaus.preds.tr,
      SPEAR.min.multi = SPEAR.min.multi.preds.tr,
      SPEAR.sd.multi = SPEAR.sd.multi.preds.tr,
      SPEAR.van.multi = SPEAR.van.multi.preds.tr,
      MOFA.multi = MOFA.multi.preds.tr,
      Lasso.multi = lasso.multi.preds.tr
    )
    colnames(pred.df.tr) <- c("SPEAR.min.gaussian",
                              "SPEAR.sd.gaussian",
                              "SPEAR.van.gaussian",
                              "MOFA.gaussian",
                              "Lasso.gaussian",
                              "SPEAR.min.multi",
                              "SPEAR.sd.multi",
                              "SPEAR.van.multi",
                              "MOFA.multi",
                              "Lasso.multi")
    pred.df.te <- data.frame(
      SPEAR.min.gaussian = SPEAR.min.gaus.preds.te,
      SPEAR.sd.gaussian = SPEAR.sd.gaus.preds.te,
      SPEAR.van.gaussian = SPEAR.van.gaus.preds.te,
      MOFA.gaussian = MOFA.gaus.preds.te,
      Lasso.gaussian = lasso.gaus.preds.te,
      SPEAR.min.multi = SPEAR.min.multi.preds.te,
      SPEAR.sd.multi = SPEAR.sd.multi.preds.te,
      SPEAR.van.multi = SPEAR.van.multi.preds.te,
      MOFA.multi = MOFA.multi.preds.te,
      Lasso.multi = lasso.multi.preds.te
    )
    colnames(pred.df.te) <- c("SPEAR.min.gaussian",
                              "SPEAR.sd.gaussian",
                              "SPEAR.van.gaussian",
                              "MOFA.gaussian",
                              "Lasso.gaussian",
                              "SPEAR.min.multi",
                              "SPEAR.sd.multi",
                              "SPEAR.van.multi",
                              "MOFA.multi",
                              "Lasso.multi")
    
    ### Misclassification Errors:
    misclass.list <- list()
    for(model in colnames(pred.df.tr)){
      # Get misclassificaiton error rate and BER:
      df <- data.frame(true = sim.data$train$Y,
                       pred = pred.df.tr[,model])
      colnames(df) <- c("true", "pred")
      cm_o <- table(df)
      cmlevels = sort(unique(df$true))
      cm<-matrix(0L,nrow=length(cmlevels),ncol=length(cmlevels))
      rownames(cm)=cmlevels
      colnames(cm)=cmlevels
      for (i in rownames(cm_o)){
        for (j in colnames(cm_o)){
          cm[i,j]=cm_o[i,j]
        }
      }
      cm_wrong<-cm
      diag(cm_wrong)<-0
      cm_wrong_sum<-apply(cm_wrong,1,sum)
      cm_sum<-apply(cm,1,sum)
      misclass_error <- round(sum(cm_wrong)/sum(cm),2)
      misclass_error_ind <- cm_wrong_sum/cm_sum
      misclass_error_ind <- na.omit(misclass_error_ind)
      ber <-round(sum(misclass_error_ind)/length(misclass_error_ind),2)
      error.tr <- misclass_error
      b.error.tr <- ber
      
      # Get misclassificaiton error rate and BER:
      df <- data.frame(true = sim.data$test$Y,
                       pred = pred.df.te[,model])
      colnames(df) <- c("true", "pred")
      cm_o <- table(df)
      cmlevels = sort(unique(df$true))
      cm<-matrix(0L,nrow=length(cmlevels),ncol=length(cmlevels))
      rownames(cm)=cmlevels
      colnames(cm)=cmlevels
      for (i in rownames(cm_o)){
        for (j in colnames(cm_o)){
          cm[i,j]=cm_o[i,j]
        }
      }
      cm_wrong<-cm
      diag(cm_wrong)<-0
      cm_wrong_sum<-apply(cm_wrong,1,sum)
      cm_sum<-apply(cm,1,sum)
      misclass_error <- round(sum(cm_wrong)/sum(cm),2)
      misclass_error_ind <- cm_wrong_sum/cm_sum
      misclass_error_ind <- na.omit(misclass_error_ind)
      ber <-round(sum(misclass_error_ind)/length(misclass_error_ind),2)
      error.te <- misclass_error
      b.error.te <- ber
      
      misclass.list[[model]] <- list(UnbalancedTrain = error.tr,
           BalancedTrain = b.error.tr,
           UnbalancedTest = error.te,
           BalancedTest = b.error.te)
    }
    
    
    ### Metadata:
    
    # CV loss:
    gaus.loss <- SPEARgaussian$get.cv.loss()
    cv.loss <- list(
      SPEARgaussian = SPEARgaussian$get.cv.loss(),
      SPEARmulti = SPEARmulti$get.cv.loss()
    )
    
    # Compile results:
    meta.list <- list(
      cv.loss = cv.loss
    )
    
    final.res[[paste0("Iteration", iteration)]] <- list(preds.tr = pred.df.tr,
                                                        preds.te = pred.df.te,
                                                        misclass.error = do.call("rbind", misclass.list),
                                                        sim.data = sim.data,
                                                        meta.list = meta.list)
}
saveRDS(final.res, paste0(save.path, "test_results_c", c.original, ".rds"))