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
  makeDiscreteParam("k", 3:9))


lrn3 = makeTuneWrapper(makeLearner("classif.kknn",predict.type="prob"), measures = auc, resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)


lrns = lrn3

outer = makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE)


bmr.hog_16_0_mosaic_1_8.cforest.lrn3 = benchmark(lrns, tasks, outer, measures = list(acc, auc,tpr,fpr), show.info = F,models = F)
