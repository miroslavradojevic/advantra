library(parallelMap)
parallelStartMulticore(20L,level="mlr.tuneParams")
library(mlr)

load(file="~/projects/ndethcmml/01_data/scspm_230_mosaic_1_8.RData")

tasks =   makeClassifTask(id = "scspm_230_mosaic_1_8", data = scspm_230_mosaic_1_8, target = "class")

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

bmr.scspm_230_mosaic_1_8.cforest.kknn = benchmark(lrns, tasks, outer, measures = list(acc, auc,tpr,fpr), show.info = T,models = F)