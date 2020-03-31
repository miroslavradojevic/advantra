library(parallelMap)
parallelStartMulticore(20L,level="mlr.tuneParams")
library(mlr)

load(file="~/projects/ndethcmml/01_data/scspm_230_mosaic_1_8.RData")

tasks =   makeClassifTask(id = "scspm_230_mosaic_1_8", data = scspm_230_mosaic_1_8, target = "class")

tasks<-removeConstantFeatures(tasks)
tasks<-normalizeFeatures(tasks, method = "standardize")


ps = makeParamSet(
  makeNumericParam("cost", lower = -12, upper = 12, trafo = function(x) 2^x)
)
ctrl = makeTuneControlGrid()
inner = makeResampleDesc("Holdout")
lrn1 = makeTuneWrapper(makeLearner("classif.svm", kernel = "linear",predict.type="prob"), measures = auc, resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)

lrns = lrn1

outer = makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE)

bmr.scspm_230_mosaic_1_8.cforest.svm = benchmark(lrns, tasks, outer, measures = list(acc, auc,tpr,fpr), show.info = T,models = F)
