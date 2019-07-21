resp <- rbinom(50, size = 1, prob = 0.5)
detpred1 <- runif(50, 0, 1)
detpred2 <- runif(50, 0, 10)
formula <- resp ~ detpred1 + detpred2
data <- data.frame(resp, detpred1, detpred2)
varmethod <- "Bootstrap"

formula <- resp ~ 1
detectionobj <- get_detection(formula = formula, data = data, 
  varmethod = varmethod)

exampledataset$detpred1 <- runif(40, 0, 1)
exampledataset$detpred2 <- runif(40, 0, 1)

exampledataset$areas <- rep(0.5, 40)
exampledataset$st <- c(rep("H", 21), rep("Low", 19))

slm_info <- slmfit(formula = counts ~ 1,
  data = exampledataset,
  xcoordcol = "xcoords", ycoordcol = "ycoords",
   estmethod = "ML",
  coordtype = "TM", detectionobj = detectionobj,
  CorModel = "Gaussian",
  areacol = "areas")

coef(slm_info)
residuals(slm_info)
summary(slm_info)
slm_info$SpatialParmEsts
slm_info$CoefficientEsts

exampledataset$st
detectionobj <- get_detection(resp ~ 1, data = data)

slm_infowithdet <- slmfit(formula = counts ~ st,
  data = exampledataset,
  xcoordcol = "xcoords", ycoordcol = "ycoords",
  estmethod = "ML", detectionobj = detectionobj,
  coordtype = "TM",
  CorModel = "Gaussian",
  areacol = NULL)
predinfo <- predict(slm_infowithdet)
FPBKoutput(predinfo, get_krigmap = TRUE)

vignettecount <- read.csv("~/Desktop/Shiny-FPBK/vignettecount.csv")


predict(slminfo)
cbind(vignettecount$Moose, vignettecount$Stratum)


coef(slm_infowithdet)
residuals(slm_infowithdet)
summary(slm_infowithdet)
slm_infowithdet$SpatialParmEsts
slm_infowithdet$CoefficientEsts

multistrat(formula = counts ~ 1,
  data = exampledataset,
  xcoordcol = "xcoords", ycoordcol = "ycoords",
  estmethod = "ML", detectionobj = NULL,
  CorModel = "Gaussian",
  areacol = NULL,
  FPBKcol = "areavar",
  stratcol = "st",
  coordtype = "TM")
?multistrat
library(FPBKPack2)
predict(object = slm_info, FPBKcol = NULL)$FPBK_Prediction
predobj <- predict(object = slm_info, FPBKcol = NULL)

tabsandstuff <- FPBKoutput(pred_info = predobj, conf_level = c(0.80, 
  0.90, 0.95),
  get_krigmap = TRUE, get_sampdetails = TRUE,
  get_variogram = TRUE,
  pointsize = 4)
  
  
get_reportdoc(output_info = tabsandstuff)

slm_info <- slmfit(formula = counts ~ pred1 + pred2,
  data = exampledataset,
  xcoordcol = "xcoords", ycoordcol = "ycoords",  coordtype = "UTM",
  estmethod = "ML", detectionobj = NULL)

predict(object = slm_info, FPBKcol = NULL,
  detinfo = c(1, 0))$FPBK_Prediction

