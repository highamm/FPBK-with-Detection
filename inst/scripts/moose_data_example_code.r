resp <- rbinom(50, size = 1, prob = 0.5)
detpred1 <- runif(50, 0, 1)
detpred2 <- runif(50, 0, 10)
formula <- resp ~ detpred1 + detpred2
data <- data.frame(resp, detpred1, detpred2)
varmethod <- "Delta"

detectionobj <- get_detection(formula = formula, data = data, 
  varmethod = varmethod)

exampledataset$detpred1 <- runif(40, 0, 1)
exampledataset$detpred2 <- runif(40, 0, 1)

slm_info <- slmfit(formula = counts ~ 1,
  data = exampledataset,
  xcoordcol = "xcoords", ycoordcol = "ycoords",  coordtype = "UTM",
   estmethod = "ML", detectionobj = detectionobj)


predict.slmfit(object = slm_info, FPBKcol = NULL) 


slm_info <- slmfit(formula = counts ~ 1,
  data = exampledataset,
  xcoordcol = "xcoords", ycoordcol = "ycoords",  coordtype = "UTM",
  estmethod = "ML", detectionobj = NULL)
predict.slmfit(object = slm_info, FPBKcol = NULL,
  detinfo = c(1, 0)) 
