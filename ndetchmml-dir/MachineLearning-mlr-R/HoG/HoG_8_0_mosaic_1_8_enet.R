library(parallelMap)
parallelStartMulticore(20L,level="mlr.tuneParams")
library(mlr)

load(file="~/projects/ndethcmml/01_data/hog_8_0_mosaic_1_8.RData")

tasks =   makeClassifTask(id = "hog_8_0_mosaic_1_8", data = hog_8_0_mosaic_1_8, target = "class")

tasks<-removeConstantFeatures(tasks)
tasks<-normalizeFeatures(tasks, method = "standardize")


ctrl = makeTuneControlGrid()
inner = makeResampleDesc("Holdout")

ps = makeParamSet(makeDiscreteParam("lambda", c(0.0001,0.001,0.01,0.1,1)),
                  makeDiscreteParam("alpha",c(0,0.15,0.25,0.35,0.5,0.65,0.75,0.85,1)))

lrn4 = makeTuneWrapper(makeLearner("classif.glmnet", predict.type = "prob"),measures = auc, resampling = inner, par.set =  ps, control = ctrl,
                       show.info = TRUE)

lrns = lrn4

outer = makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE)

bmr.hog_8_0_mosaic_1_8.cforest.enet = benchmark(lrns, tasks, outer, measures = list(acc, auc,tpr,fpr), show.info = T,models = F)