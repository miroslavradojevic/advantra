library(parallelMap)
parallelStartMulticore(20L,level="mlr.tuneParams")
library(mlr)
library(foreign)

data_rnd_CHARM = read.arff(file='~/ndethcmml/rnd_CHARM_D.S.OPOS_500_1.0_0.5.arff')
names.data<-paste('V',c(1:1059),sep='')
colnames(data_rnd_CHARM)<-make.names(c(names.data,'class'))

data_rnd_CHARM_SIFT_230 = read.arff(file='~/ndethcmml/rnd_CHARM-SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_230_11.arff')
names.data<-paste('V',c(1:1289),sep='')
colnames(data_rnd_CHARM_SIFT_230)<-make.names(c(names.data,'class'))

data_trainfeats_230 = read.arff(file='~/ndethcmml/rnd_SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_230_11.arff')
names.data<-paste('V',c(1:230),sep='')
colnames(data_trainfeats_230)<-make.names(c(names.data,'class'))


data_rnd_CHARM<-normalizeFeatures(data_rnd_CHARM, method = "standardize")
data_rnd_CHARM_SIFT_230<-normalizeFeatures(data_rnd_CHARM_SIFT_230, method = "standardize")
data_trainfeats_230<-normalizeFeatures(data_trainfeats_230, method = "standardize")

data_rnd_CHARM<-removeConstantFeatures(data_rnd_CHARM)
data_rnd_CHARM_SIFT_230<-removeConstantFeatures(data_rnd_CHARM_SIFT_230)
data_trainfeats_230<-removeConstantFeatures(data_trainfeats_230)

tasks = list(
  makeClassifTask(id = "CHARM", data = data_rnd_CHARM, target = "class"),
  makeClassifTask(id = "CHARM_SIFT_230", data = data_rnd_CHARM_SIFT_230, target = "class"),
  makeClassifTask(id = "trainfeats_230", data = data_trainfeats_230, target = "class")
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

ctrl = makeTuneControlGrid()
inner = makeResampleDesc("Holdout")
lrn2 = makeTuneWrapper(makeLearner("classif.randomForest",predict.type="prob",config = list(importance=TRUE,ntree=1000L)), measures = auc, resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)
ps = makeParamSet(
  makeDiscreteParam("k", 3:9))

ctrl = makeTuneControlGrid()
inner = makeResampleDesc("Holdout")
lrn3 = makeTuneWrapper(makeLearner("classif.kknn",predict.type="prob"), measures = auc, resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)

ps = makeParamSet(makeDiscreteParam("lambda", c(0.0001,0.001,0.01,0.1,1)),
                  makeDiscreteParam("alpha",c(0,0.15,0.25,0.35,0.5,0.65,0.75,0.85,1)))
ctrl = makeTuneControlGrid()
inner = makeResampleDesc("Holdout")
lrn4 = makeTuneWrapper(makeLearner("classif.glmnet", predict.type = "prob"),measures = auc, resampling = inner, par.set =  ps, control = ctrl,
                       show.info = TRUE)


lrns = list(lrn1,lrn2,lrn3,lrn4)

outer = list(makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE))


res = benchmark(lrns, tasks, outer, measures = list(acc, auc,tpr,fpr), show.info = TRUE,models = TRUE)
