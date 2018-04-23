library(parallelMap)
parallelStartMulticore(20L,level="mlr.tuneParams")
library(mlr)
library(foreign)


data_rnd_CHARM_SIFT_230 = read.arff(file='~/projects/ndethcmml/01_data/rnd_CHARM-SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_230_11.arff')
names.data<-paste('V',c(1:1289),sep='')
colnames(data_rnd_CHARM_SIFT_230)<-make.names(c(names.data,'class'))

data_rnd_CHARM_SIFT_200 = read.arff(file='~/projects/ndethcmml/01_data/rnd_CHARM-SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_200_11.arff')
names.data<-paste('V',c(1:1259),sep='')
colnames(data_rnd_CHARM_SIFT_200)<-make.names(c(names.data,'class'))

data_rnd_CHARM_SIFT_40 = read.arff(file='~/projects/ndethcmml/01_data/rnd_CHARM-SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_40_11.arff')
names.data<-paste('V',c(1:1099),sep='')
colnames(data_rnd_CHARM_SIFT_40)<-make.names(c(names.data,'class'))

data_rnd_CHARM_SIFT_20 = read.arff(file='~/projects/ndethcmml/01_data/rnd_CHARM-SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_20_11.arff')
names.data<-paste('V',c(1:1079),sep='')
colnames(data_rnd_CHARM_SIFT_20)<-make.names(c(names.data,'class'))

data_rnd_CHARM_SIFT_150 = read.arff(file='~/projects/ndethcmml/01_data/rnd_CHARM-SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_150_11.arff')
names.data<-paste('V',c(1:1209),sep='')
colnames(data_rnd_CHARM_SIFT_150)<-make.names(c(names.data,'class'))
cat('cargando datos\n')
data_rnd_CHARM_SIFT_80 = read.arff(file='~/projects/ndethcmml/01_data/rnd_CHARM-SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_80_11.arff')
names.data<-paste('V',c(1:1139),sep='')
colnames(data_rnd_CHARM_SIFT_80)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

data_rnd_CHARM_SIFT_60 = read.arff(file='~/projects/ndethcmml/01_data/rnd_CHARM-SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_60_11.arff')
names.data<-paste('V',c(1:1119),sep='')
colnames(data_rnd_CHARM_SIFT_60)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

data_rnd_CHARM_SIFT_100 = read.arff(file='~/projects/ndethcmml/01_data/rnd_CHARM-SIFT.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_100_11.arff')
names.data<-paste('V',c(1:1159),sep='')
colnames(data_rnd_CHARM_SIFT_100)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

data_trainfeats_100 = read.arff(file='~/projects/ndethcmml/01_data/rnd_trainfeats.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_100_11.arff')
names.data<-paste('V',c(1:100),sep='')
colnames(data_trainfeats_100)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

data_trainfeats_60 = read.arff(file='~/projects/ndethcmml/01_data/rnd_trainfeats.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_60_11.arff')
names.data<-paste('V',c(1:60),sep='')
colnames(data_trainfeats_60)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

data_trainfeats_80 = read.arff(file='~/projects/ndethcmml/01_data/rnd_trainfeats.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_80_11.arff')
names.data<-paste('V',c(1:80),sep='')
colnames(data_trainfeats_80)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

data_trainfeats_150 = read.arff(file='~/projects/ndethcmml/01_data/rnd_trainfeats.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_150_11.arff')
names.data<-paste('V',c(1:150),sep='')
colnames(data_trainfeats_150)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

data_trainfeats_20 = read.arff(file='~/projects/ndethcmml/01_data/rnd_trainfeats.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_20_11.arff')
names.data<-paste('V',c(1:20),sep='')
colnames(data_trainfeats_20)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

data_trainfeats_40 = read.arff(file='~/projects/ndethcmml/01_data/rnd_trainfeats.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_40_11.arff')
names.data<-paste('V',c(1:40),sep='')
colnames(data_trainfeats_40)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

data_trainfeats_200 = read.arff(file='~/projects/ndethcmml/01_data/rnd_trainfeats.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_200_11.arff')
names.data<-paste('V',c(1:200),sep='')
colnames(data_trainfeats_200)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

data_trainfeats_230 = read.arff(file='~/projects/ndethcmml/01_data/rnd_trainfeats.comb.D.S.OPOS.NW.NPOS_comb1_500_1.0_0.5_230_11.arff')
names.data<-paste('V',c(1:230),sep='')
colnames(data_trainfeats_230)<-make.names(c(names.data,'class'))
cat('cargando datos\n')

tasks = list(
makeClassifTask(id = "CHARM_SIFT_230", data = data_rnd_CHARM_SIFT_230, target = "class"),
makeClassifTask(id = "CHARM_SIFT_200", data = data_rnd_CHARM_SIFT_200, target = "class"),
makeClassifTask(id = "CHARM_SIFT_40", data = data_rnd_CHARM_SIFT_40, target = "class"),
makeClassifTask(id = "CHARM_SIFT_20", data = data_rnd_CHARM_SIFT_20, target = "class"),
  makeClassifTask(id = "CHARM_SIFT_150", data = data_rnd_CHARM_SIFT_150, target = "class"),
  makeClassifTask(id = "CHARM_SIFT_80", data = data_rnd_CHARM_SIFT_80, target = "class"),
  makeClassifTask(id = "CHARM_SIFT_60", data = data_rnd_CHARM_SIFT_60, target = "class"),
  makeClassifTask(id = "CHARM_SIFT_100", data = data_rnd_CHARM_SIFT_100, target = "class"),
  makeClassifTask(id = "trainfeats_100", data = data_trainfeats_100, target = "class"),
  makeClassifTask(id = "trainfeats_20", data = data_trainfeats_20, target = "class"),
  makeClassifTask(id = "trainfeats_40", data = data_trainfeats_40, target = "class"),
  makeClassifTask(id = "trainfeats_200", data = data_trainfeats_200, target = "class"),
  makeClassifTask(id = "trainfeats_230", data = data_trainfeats_230, target = "class"),
  makeClassifTask(id = "trainfeats_60", data = data_trainfeats_60, target = "class"),
  makeClassifTask(id = "trainfeats_80", data = data_trainfeats_80, target = "class"),
  makeClassifTask(id = "trainfeats_150", data = data_trainfeats_150, target = "class")
)

ps = makeParamSet(
  makeNumericParam("C", lower = -12, upper = 12, trafo = function(x) 2^x),
  makeNumericParam("sigma", lower = -12, upper = 12, trafo = function(x) 2^x)
)
ctrl = makeTuneControlGrid()
inner = makeResampleDesc("Holdout")
lrn1 = makeTuneWrapper(makeLearner("classif.ksvm",predict.type="prob"), resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)

ps = makeParamSet(
  makeDiscreteParam("mtry", values=c(5,10,15,20,50,100)),
  makeDiscreteParam("ntree", values=c(50,100,500)))

ctrl = makeTuneControlGrid()
lrn2 = makeTuneWrapper(makeLearner("classif.randomForest",predict.type="prob"), resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)
ps = makeParamSet(
  makeDiscreteParam("k", 3:5))

ctrl = makeTuneControlGrid()
lrn3 = makeTuneWrapper(makeLearner("classif.kknn",predict.type="prob"), resampling = inner, par.set = ps, control = ctrl,
                       show.info = TRUE)

ps = makeParamSet(makeDiscreteParam("lambda", c(0.001,0.01,0.1,1)),
                  makeDiscreteParam("alpha",c(0,0.25,0.5,0.75,1)))
ctrl = makeTuneControlGrid()
lrn4 = makeTuneWrapper(makeLearner("classif.glmnet", predict.type = "prob"), resampling = inner, par.set =  ps, control = ctrl,
                       show.info = TRUE)


lrns = list(lrn1,lrn2,lrn3,lrn4)

outer = list(makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE),
             makeResampleDesc("RepCV",reps=3,folds=10,stratify=TRUE))
 
res = benchmark(lrns, tasks, outer, measures = list(acc, auc), show.info = TRUE)


