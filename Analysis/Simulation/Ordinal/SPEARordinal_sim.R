library(SPEAR)
library(tidyverse)

# Path to save results to:
save.path <- "~/SPEAR_ordinal_results"

final.results <- list()

for(rep in 1:10){
  total.results <- list()
  for(idx in 1:3){
    
      sim.data <- SPEAR::simulate_data(
        family = "ordinal",
        N = 500,
        P = 200,
        ordinal.centers = c(-1.5, -1.3, -1, 0, 1, 1.3, 1.5),
        c1 = c(.1, .5, 1, 3)[idx],
        c2 = 10,
        seed = 122 + rep
      )
      
      SPEARgaussian <- SPEAR::make.SPEARobject(X = sim.data$train$Xlist,
                                               Y = sim.data$train$Y,
                                               family = "gaussian")
      SPEARgaussian$add.data(X = sim.data$test$Xlist, Y = sim.data$test$Y, name = "test")
      SPEARgaussian$run.cv.spear()
      SPEARgaussian$cv.evaluate()
      SPEARgaussian$set.weights(method = "sd")
      
      
      SPEARordinal <- SPEAR::make.SPEARobject(X = sim.data$train$Xlist,
                                              Y = sim.data$train$Y,
                                              family = "ordinal")
      SPEARordinal$add.data(X = sim.data$test$Xlist, Y = sim.data$test$Y, name = "test")
      SPEARordinal$run.cv.spear()
      SPEARordinal$cv.evaluate()
      SPEARordinal$set.weights(method = "sd")
      
      
      SPEARmulti <- SPEAR::make.SPEARobject(X = sim.data$train$Xlist,
                                            Y = sim.data$train$Y.onehot,
                                            family = "multinomial")
      SPEARmulti$add.data(X = sim.data$test$Xlist, Y = sim.data$test$Y.onehot, name = "test")
      SPEARmulti$run.cv.spear()
      SPEARmulti$cv.evaluate()
      SPEARmulti$set.weights(method = "sd")
      
      
      errors <- list()
      
      SPEAR.gaus.preds.te <- SPEARgaussian$get.predictions(data = "test")  + SPEARgaussian$fit$intercepts.scaled[[1]][SPEARgaussian$options$current.weight.idx,]
      SPEAR.gaus.preds.te <- round(SPEAR.gaus.preds.te, 0)
      SPEAR.gaus.preds.te[SPEAR.gaus.preds.te < 0] <- 0
      SPEAR.gaus.preds.te[SPEAR.gaus.preds.te > 6] <- 6
      
      df <- data.frame(true.vals = sim.data$test$Y,
                       pred.vals = SPEAR.gaus.preds.te)
      colnames(df) <- c("true", "pred")
      
      cm_o <- table(df)
      cmlevels = 0:6
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
      
      errors$gaussian <- list(misclass_error = misclass_error, ber = ber)
      o.errors <- SPEARordinal$get.misclassification(data = "test")
      errors$ordinal <- list(misclass_error = o.errors$misclass_error, ber = o.errors$balanced_misclass_error)
      m.errors <- SPEARmulti$get.misclassification(data = "test")
      errors$multi <- list(misclass_error = m.errors$misclass_error, ber = m.errors$balanced_misclass_error)
      
      # Gaussian:
      true <- data.frame(true = sim.data$test$Y)
      preds <- data.frame(
        SPEARgaussian = SPEAR.gaus.preds.te,
        SPEARordinal = SPEARordinal$get.predictions(data = "test")[[1]]$class.predictions[,1],
        SPEARmulti = SPEARmulti$get.predictions(data = "test")$class.predictions[,1],
        TrueVals = sim.data$test$Y,
        true = sim.data$test$Y
      )
      colnames(preds)[1] <- "SPEARgaussian"
      
      total.results[[c("low", "med", "high")[idx]]] <- list(
        errors = errors,
        df = preds
      )
  }
  final.results[[rep]] <- total.results
}

saveRDS(final.results, paste0(save.path, "/finalresults.rds"))

