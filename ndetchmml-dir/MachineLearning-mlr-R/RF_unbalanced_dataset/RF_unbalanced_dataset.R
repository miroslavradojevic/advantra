library(parallelMap)
parallelStartMulticore(20L,level="mlr.tuneParams")
library(mlr)

library(foreign)

charm.original = read.arff(file='~/projects/ndethcmml/01_data/CHARM_unbalanced.arff')

names.data<-paste('V',c(1:1059),sep='')
colnames(charm.original)<-make.names(c(names.data,'class'))

charm.sift.original = read.arff(file='~/projects/ndethcmml/01_data/CHARM-SIFT_unbalanced.arff')

names.data<-paste('V',c(1:1289),sep='')
colnames(charm.sift.original)<-make.names(c(names.data,'class'))

sift.original = read.arff(file='~/projects/ndethcmml/01_data/SIFT_unbalanced.arff')

names.data<-paste('V',c(1:230),sep='')
colnames(sift.original)<-make.names(c(names.data,'class'))


charm.task=makeClassifTask(id = "CHARM_unbal", data = charm.original, target = "class")
charm.sift.task=makeClassifTask(id = "CHARM_SIFT_unbal", data = charm.sift.original, target = "class")
sift.task=makeClassifTask(id = "SIFT_unbal", data = sift.original, target = "class")

charm.task=removeConstantFeatures(charm.task)
charm.sift.task=removeConstantFeatures(charm.sift.task)
sift.task=removeConstantFeatures(sift.task)

charm.task<-normalizeFeatures(charm.task, method = "standardize")
charm.sift.task<-normalizeFeatures(charm.sift.task, method = "standardize")
sift.task<-normalizeFeatures(sift.task, method = "standardize")


tasks = list(
  charm.task,
  charm.sift.task,
  sift.task
)

ctrl = makeTuneControlGrid()
inner = makeResampleDesc("Holdout")

ps = makeParamSet(
  makeDiscreteParam("nodesize", values=c(1:5)),
  makeDiscreteParam("mtry", values=c(5:36)))

lrn2 = makeTuneWrapper(makeLearner("classif.randomForest",predict.type="prob",par.vals = list(importance=TRUE,ntree=1000L)), measures = auc, resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)
lrns = lrn2

outer = list(
			makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
			makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
			makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE)
			)

bmr.rf.unbalance = benchmark(lrns, tasks, outer, measures = list(acc, auc,tpr,fpr), show.info = F,models = F)

