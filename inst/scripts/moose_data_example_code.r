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

exampledataset$areas <- rep(0.5, 40)
exampledataset$st <- c(rep("H", 25), rep("Low", 15))

slm_info <- slmfit(formula = counts ~ 1,
  data = exampledataset,
  xcoordcol = "xcoords", ycoordcol = "ycoords", coordtype = "UTM",
   estmethod = "ML", detectionobj = NULL,
  CorModel = "Gaussian",
  areacol = "areas")
slm_info$SpatialParmEsts
slm_info$CoefficientEsts

predict(object = slm_info, FPBKcol = NULL)$FPBK_Prediction
predobj <- predict(object = slm_info, FPBKcol = NULL)

FPBKoutput(pred_info = predobj, conf_level = c(0.80, 
  0.90, 0.95),
  get_krigmap = TRUE, get_sampdetails = TRUE,
  get_variogram = TRUE, get_report = TRUE,
  pointsize = 4)
  
  
slm_info <- slmfit(formula = counts ~ pred1 + pred2,
  data = exampledataset,
  xcoordcol = "xcoords", ycoordcol = "ycoords",  coordtype = "UTM",
  estmethod = "ML", detectionobj = NULL)

predict(object = slm_info, FPBKcol = NULL,
  detinfo = c(1, 0))$FPBK_Prediction

