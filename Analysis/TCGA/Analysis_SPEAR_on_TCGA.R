### Load Libraries:
# lasso
library(glmnet)
library(matrixStats)
# MOFA
require(parallel)
library(reticulate)
library(MOFA2)
predict <- stats::predict
# DIABLO
library(mixOmics)
# SPEAR
library(SPEAR)

###########################################
#     REPRODUCIBILITY NOTE FOR USER:      #
#     CHANGE THE PYTHON PATH HERE:        #
python_path <- "/opt/anaconda3/bin/python"#
###########################################
use_python(python_path, required = TRUE)

# Path to save results to:
save.path <- "~/SPEAR_TCGA/results/"
tcga.data.location <- "~/SPEAR_TCGA/data/"

# initial seed:
seed <- 124

# Load data:
load(paste0(tcga.data.location, "Singh_TCGA_TCGA.normalised.mixDIABLO.RData")) 
load(paste0(tcga.data.location, "Singh_TCGA_trainTestDatasetsNormalized.RDATA"))

# make into lists
full_train<-list(meth=methTrain0,mirna=mirnaTrain0,mrna=mrnaTrain0,pam50=pam50Train0,clin=clinTrain0)
full_test<-list(meth=methTest0,mirna=mirnaTest0,mrna=mrnaTest0,pam50=pam50Test0,clin=clinTest0)

### Feature filtering
mrna_full_var<-data.frame(row.names=colnames(full_train$mrna),colVars(full_train$mrna))
mrna_full_mean<-data.frame(row.names=colnames(full_train$mrna),colMeans(full_train$mrna))
mrna_full_var$hold<-1
colnames(mrna_full_var)<-c('Var','hold')
mrna_full_mean$hold<-1
colnames(mrna_full_mean)<-c('Mean','hold')
mrna_full_var$Mean<-mrna_full_mean$Mean
mrna_full_var<-mrna_full_var[order(mrna_full_var$Var,decreasing=TRUE),]

quant_mrna<-quantile(mrna_full_var$Var,c(0.2))
mrna_keep<-mrna_full_var[mrna_full_var$Var>=quant_mrna[1],]

meth_full_var<-data.frame(row.names=colnames(full_train$meth),colVars(full_train$meth))
meth_full_var$hold<-1
colnames(meth_full_var)<-c('Var','hold')
meth_full_var<-meth_full_var[order(meth_full_var$Var,decreasing=TRUE),]
quant_meth<-quantile(meth_full_var$Var,c(0.2))
meth_keep<-meth_full_var[meth_full_var$Var>=quant_meth[1],]


mirna_full_var<-data.frame(row.names=colnames(full_train$mirna),colVars(full_train$mirna))
mirna_full_var$hold<-1
colnames(mirna_full_var)<-c('Var','hold')
mirna_full_var<-mirna_full_var[order(mirna_full_var$Var,decreasing=TRUE),]
quant_mirna<-quantile(mirna_full_var$Var,c(0.2))
mirna_keep<-mirna_full_var[mirna_full_var$Var>=quant_mirna[1],]

full_train_f<-list(mrna=mrnaTrain0[,rownames(mrna_keep)],
                   mirna=mirnaTrain0[,rownames(mirna_keep)],
                   meth=methTrain0[,rownames(meth_keep)],
                   clin=clinTrain0, subtype=data.train$subtype)
lapply(full_train_f, dim)

full_test_f<-list(mrna=mrnaTest0[,rownames(mrna_keep)],
                  mirna=mirnaTest0[,rownames(mirna_keep)],
                  meth=methTest0[,rownames(meth_keep)],
                  clin=clinTest0, subtype=data.test$subtype)
lapply(full_train_f, dim)

# scale
full_train_f$mrna<-scale(full_train_f$mrna)
full_train_f$mirna<-scale(full_train_f$mirna)
full_train_f$meth<-scale(full_train_f$meth)

full_test_f$mrna<-scale(full_test_f$mrna)
full_test_f$mirna<-scale(full_test_f$mirna)
full_test_f$meth<-scale(full_test_f$meth)

### Prep for multinomial subtype
full_train_f$subtype<-factor(full_train_f$subtype, levels=c('Basal','Her2','LumB','LumA'))
full_train_f$num<-sapply(full_train_f$subtype,unclass)-1

full_test_f$subtype<-factor(full_test_f$subtype, levels=c('Basal','Her2','LumB','LumA'))
full_test_f$num<-sapply(full_test_f$subtype,unclass)-1

Y_mult_train <-data.frame(matrix(0,379,4))
colnames(Y_mult_train)<-c('Basal','Her2','LumB','LumA')

r_basal<-which(full_train_f$subtype=='Basal')
r_HER2<-which(full_train_f$subtype=='Her2')
r_LumB<-which(full_train_f$subtype=='LumB')
r_LumA<-which(full_train_f$subtype=='LumA')

Y_mult_train[r_basal,'Basal'] =1
Y_mult_train[r_HER2,'Her2']=1
Y_mult_train[r_LumB,'LumB']=1
Y_mult_train[r_LumA,'LumA']=1

rownames(Y_mult_train)<-rownames(full_train_f$mrna)

Y_mult_test <-data.frame(matrix(0,610,4))
colnames(Y_mult_test)<-c('Basal','Her2','LumB','LumA')

r_basal<-which(full_test_f$subtype=='Basal')
r_HER2<-which(full_test_f$subtype=='Her2')
r_LumB<-which(full_test_f$subtype=='LumB')
r_LumA<-which(full_test_f$subtype=='LumA')

Y_mult_test[r_basal,'Basal'] =1
Y_mult_test[r_HER2,'Her2']=1
Y_mult_test[r_LumB,'LumB']=1
Y_mult_test[r_LumA,'LumA']=1

rownames(Y_mult_test)<-rownames(full_test_f$mrna)

full_train_f$subtype_mult<-Y_mult_train
full_test_f$subtype_mult<-Y_mult_test

N = nrow(Y_mult_train)
nfolds = 10
foldid = sample(rep(1:nfolds, ceiling(N/nfolds)), N, replace = F)

X <- list(mrna = full_train_f$mrna,
          mirna = full_train_f$mirna,
          meth = full_train_f$meth)
Y <- matrix(full_train_f$subtype, ncol = 1)
Y.multi <- full_train_f$subtype_mult
X.te <- list(mrna = full_test_f$mrna,
          mirna = full_test_f$mirna,
          meth = full_test_f$meth)
Y.te <- matrix(full_test_f$subtype, ncol = 1)
Y.te.multi <- full_test_f$subtype_mult
get_confusion_matrix_from_predictions <- function(df){
  # NOTE: first column must be true, second column must be predicted
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
  cm <- as.table(cm)
  names(dimnames(cm)) <- c("True", "Pred")
  return(cm)
}
get_errors <- function(pred.df){
  errors <- list()
  # unbalanced misclassification error:
  conf.mat <- get_confusion_matrix_from_predictions(pred.df)
  cm_wrong<-conf.mat
  diag(cm_wrong)<-0
  cm_wrong_sum<-apply(cm_wrong,1,sum)
  cm_sum<-apply(conf.mat,1,sum)
  misclass_error <- round(sum(cm_wrong)/sum(conf.mat),2)
  misclass_error_ind <- cm_wrong_sum/cm_sum
  misclass_error_ind <- na.omit(misclass_error_ind)
  
  # balanced misclassification error:
  ber <-round(sum(misclass_error_ind)/length(misclass_error_ind),2)
  
  # individual misclassification errors (per class):
  misclass_error_ind <- round(misclass_error_ind,2)
  
  errors$misclass_error <- misclass_error
  errors$bal_misclass_error <- ber
  errors$misclass_error_ind <- misclass_error_ind
  
  return(errors)
}



  
# Get arguments:
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
  # Data parameters:
  model.choice = "MOFA"
  num.factors.in.model = fact.loop
} else {
  # Data parameters:
  model.choice = args[1] # Name of model to test
  num.factors.in.model = as.numeric(args[2]) # if applicable, how many factors?
}

if(model.choice == "SPEAR"){
  
  SPEARmulti <- make.SPEARobject(X = X,
                                 Y = Y.multi,
                                 family = "multinomial",
                                 num_factors = num.factors.in.model,
                                 seed = seed)
  SPEARmulti$add.data(X = X.te, Y = Y.te.multi, name = "test")
  SPEARmulti$run.cv.spear(fold.ids = foldid, parallel.method = "parLapply")
  SPEARmulti$cv.evaluate()
  SPEARmulti$save.model(file = paste0("SPEAR_full_", num.factors.in.model, ".rds"))
  
  # SD:
  SPEARmulti$set.weights(method = "sd")
  preds.df <- data.frame(true = apply(Y.te.multi, 1, which.max),
                         pred = SPEARmulti$get.predictions(data = "test")$class.predictions[,1] + 1)
  colnames(preds.df) <- c("true", "pred")
  errors <- get_errors(preds.df)
  saveRDS(errors, paste0(save.path, "SPEARsd_", num.factors.in.model, ".rds"))
  
  # min:
  SPEARmulti$set.weights(method = "min")
  preds.df <- data.frame(true = apply(Y.te.multi, 1, which.max),
                         pred = SPEARmulti$get.predictions(data = "test")$class.predictions[,1] + 1)
  colnames(preds.df) <- c("true", "pred")
  errors <- get_errors(preds.df)
  saveRDS(errors, paste0(save.path, "SPEARmin_", num.factors.in.model, ".rds"))
  
  
} else if(model.choice == "Lasso"){
    lasso.multi_fit = cv.glmnet(x = do.call("cbind", X), y = as.matrix(Y.multi), foldid = foldid, family = "multinomial")
    # Predict:
    lasso.preds.te = predict(lasso.multi_fit, do.call("cbind", X.te), s = "lambda.min")
    preds.df <- data.frame(true = apply(Y.te.multi, 1, which.max),
                           pred = apply(lasso.preds.te, 1, which.max))
    colnames(preds.df) <- c("true", "pred")
    errors <- get_errors(preds.df)
    saveRDS(errors, paste0(save.path, "Lasso_4.rds"))
    
} else if(model.choice == "MOFA"){
    # Run MOFA:
    train.data <- list()
    train.data$X <- X
    test.data <- list()
    test.data$X <- X.te
    x.mofa <- list()
    for(d in 1:length(train.data$X)){
      x.mofa[[d]] <- t(train.data$X[[d]])
    }
    MOFAobject <- MOFA2::create_mofa(x.mofa)
    # MOFA+ specific parameters:
    data_opts <- get_default_data_options(MOFAobject)
    model_opts <- get_default_model_options(MOFAobject)
    model_opts$num_factors <- num.factors.in.model
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
    cat("DONE WITH MOFA\n")
    if(ncol(MOFA_res@expectations$Z[[1]]) < 2){
      MOFA.factors <- cbind(MOFA_res@expectations$Z[[1]], MOFA_res@expectations$Z[[1]])
      colnames(MOFA.factors) <- c("Factor1", "Factor2")
    } else {
      MOFA.factors <- MOFA_res@expectations$Z[[1]]
    }
    MOFA.coefs = matrix(NA, ncol = ncol(MOFA.factors), nrow = ncol(do.call("cbind", train.data$X)))
    for(k in 1:ncol(MOFA.factors)){
      cat(k, ":", ncol(MOFA.factors), "...\n")
      tmp =  glmnet::cv.glmnet(x = do.call("cbind", train.data$X), y = MOFA.factors[,k])
      MOFA.coefs[,k] = coef(tmp, s = "lambda.min")[-1]
    }
    MOFA.factors.tr = do.call("cbind", train.data$X) %*% MOFA.coefs
    MOFA.factors.te = do.call("cbind", test.data$X) %*% MOFA.coefs
    
    # Predict:
    MOFA_lasso_res = cv.glmnet(MOFA.factors.tr, Y, family = "multinomial", foldid = foldid)
    MOFA.preds.te = predict(MOFA_lasso_res, MOFA.factors.te, s = "lambda.min")
    
    # Get misclassificaiton error rate and BER:
    preds.df <- data.frame(true = apply(Y.te.multi, 1, which.max),
                     pred = apply(MOFA.preds.te[,,1], 1, which.max))
    colnames(preds.df) <- c("true", "pred")
    # for some reason, mofa switches 3+4...
    preds.df$pred[preds.df$pred == 3] <- 5
    preds.df$pred[preds.df$pred == 4] <- 3
    preds.df$pred[preds.df$pred == 5] <- 4
    
    errors <- get_errors(preds.df)
    saveRDS(errors, paste0(save.path, "MOFA_", num.factors.in.model, ".rds"))
} else if(model.choice == "DIABLO"){
  
  ncomp = 3
  # Taken from DIABLO supplementary information
  list.keepX <- list(mrna = c(20, 5, 20),
                     mirna = c(20, 5, 20),
                     meth = c(15, 5, 5))
  
  design = matrix(0.1, nrow = 4, ncol = 4)
  design[,ncol(design)] <- 1
  design[nrow(design),] <- 1
  diag(design) = 0
  diablo.res = block.splsda(X = X, Y = factor(Y), ncomp = ncomp, 
                            keepX = list.keepX, design = design)
  
  # Testing:
  predict.diablo.tr = predict(diablo.res, newdata = X)
  predict.diablo.te = predict(diablo.res, newdata = X.te)
  
  # Probabilities:
  predict.diablo.probs.tr <- predict.diablo.tr$WeightedPredict[,,ncomp]
  predict.diablo.probs.te <- predict.diablo.te$WeightedPredict[,,ncomp]
  #norm.tr <- t(apply(predict.diablo.probs.tr, 1, function(row){(row - min(row))/(max(row) - min(row))/sum((row - min(row))/(max(row) - min(row)))}))
  preds.df <- data.frame(true = Y.te,
                   pred = predict.diablo.te$WeightedVote$mahalanobis.dist[,ncomp])
  colnames(preds.df) <- c("true", "pred")
  errors <- get_errors(preds.df)
  saveRDS(errors, paste0(save.path, "DIABLO_4.rds"))
}
