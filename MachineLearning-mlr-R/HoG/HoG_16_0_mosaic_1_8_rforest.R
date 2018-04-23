library(parallelMap)
parallelStartMulticore(20L,level="mlr.tuneParams")
library(mlr)

load(file="~/projects/ndethcmml/01_data/hog_16_0_mosaic_1_8.RData")

tasks =   makeClassifTask(id = "hog_16_0_mosaic_1_8", data = hog_16_0_mosaic_1_8, target = "class")

tasks<-removeConstantFeatures(tasks)
tasks<-normalizeFeatures(tasks, method = "standardize")


ctrl = makeTuneControlGrid()
inner = makeResampleDesc("Holdout")


ps = makeParamSet(
  makeDiscreteParam("nodesize", values=c(1:5)),
  makeDiscreteParam("mtry", values=c(5:36)))


lrn2 = makeTuneWrapper(makeLearner("classif.randomForest",predict.type="prob",par.vals = list(importance=TRUE,ntree=1000L)), measures = auc, resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)
lrns = lrn2

outer = makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE)


bmr.hog_16_0_mosaic_1_8.cforest.rf = benchmark(lrns, tasks, outer, measures = list(acc, auc,tpr,fpr), show.info = F,models = F)
