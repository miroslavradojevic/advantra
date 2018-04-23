library(parallelMap)
parallelStartMulticore(20L,level="mlr.tuneParams")
library(mlr)

load(file="~/projects/ndethcmml/01_data/hog_4_0_mosaic_1_8.RData")

tasks =   makeClassifTask(id = "hog_4_0_mosaic_1_8", data = hog_4_0_mosaic_1_8, target = "class")

tasks<-removeConstantFeatures(tasks)
tasks<-normalizeFeatures(tasks, method = "standardize")


ps = makeParamSet(
  makeNumericParam("cost", lower = -12, upper = 12, trafo = function(x) 2^x),
  makeNumericParam("gamma", lower = -12, upper = 12, trafo = function(x) 2^x)
)
ctrl = makeTuneControlGrid()
inner = makeResampleDesc("Holdout")
lrn1 = makeTuneWrapper(makeLearner("classif.svm",predict.type="prob"), measures = auc, resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)

ps = makeParamSet(
  makeDiscreteParam("nodesize", values=c(1:5)),
  makeDiscreteParam("mtry", values=c(5:36)))

lrn2 = makeTuneWrapper(makeLearner("classif.randomForest",predict.type="prob",par.vals = list(importance=TRUE,ntree=1000L)), measures = auc, resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)
ps = makeParamSet(
  makeDiscreteParam("k", 3:9))

lrn3 = makeTuneWrapper(makeLearner("classif.kknn",predict.type="prob"), measures = auc, resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)

ps = makeParamSet(makeDiscreteParam("lambda", c(0.0001,0.001,0.01,0.1,1)),
                  makeDiscreteParam("alpha",c(0,0.15,0.25,0.35,0.5,0.65,0.75,0.85,1)))

lrn4 = makeTuneWrapper(makeLearner("classif.glmnet", predict.type = "prob"),measures = auc, resampling = inner, par.set =  ps, control = ctrl,
                       show.info = TRUE)

lrns = list(lrn1,lrn2,lrn3,lrn4)

outer = makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE)

bmr.hog_4_0_mosaic_1_8.cforest = benchmark(lrns, tasks, outer, measures = list(acc, auc,tpr,fpr), show.info = F,models = F)