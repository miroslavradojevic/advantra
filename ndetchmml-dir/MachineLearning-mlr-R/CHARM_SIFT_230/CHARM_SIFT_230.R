library(parallelMap)
parallelStartMulticore(20L,level="mlr.tuneParams")
library(mlr)
library(foreign)

data_rnd_CHARM_SIFT_230 = read.arff(file='~/projects/ndethcmml/01_data/rnd_CHARM-SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_230_11.arff')

names.data<-paste('V',c(1:1289),sep='')
colnames(data_rnd_CHARM_SIFT_230)<-make.names(c(names.data,'class'))

tasks =   makeClassifTask(id = "CHARM_SIFT_230", data = data_rnd_CHARM_SIFT_230, target = "class")

tasks<-removeConstantFeatures(tasks)

filt_ttest = generateFilterValuesData(tasks, method = "cforest.importance")

filtered.task.25 = filterFeatures(tasks,fval=filt_ttest,abs=25)
filtered.task.100 = filterFeatures(tasks,fval=filt_ttest,abs=100)
filtered.task.200 = filterFeatures(tasks,fval=filt_ttest,abs=200)
filtered.task.600 = filterFeatures(tasks,fval=filt_ttest,abs=600)

filtered.task.25.norm<-normalizeFeatures(filtered.task.25, method = "standardize")
filtered.task.100.norm<-normalizeFeatures(filtered.task.100, method = "standardize")
filtered.task.200.norm<-normalizeFeatures(filtered.task.200, method = "standardize")
filtered.task.600.norm<-normalizeFeatures(filtered.task.600, method = "standardize")

filtered.task.25.norm$task.desc$id<-'CHARM_SIFT_230_25'
filtered.task.100.norm$task.desc$id<-'CHARM_SIFT_230_100'
filtered.task.200.norm$task.desc$id<-'CHARM_SIFT_230_200'
filtered.task.600.norm$task.desc$id<-'CHARM_SIFT_230_600'

tasks = list(
  filtered.task.25.norm,
  filtered.task.100.norm,
  filtered.task.200.norm,
  filtered.task.600.norm
)

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

outer = list(makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE))


bmr.pub.cforest = benchmark(lrns, tasks, outer, measures = list(acc, auc,tpr,fpr), show.info = TRUE,models = TRUE)
