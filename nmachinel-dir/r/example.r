## clear 
rm(list = ls());

this.dir <- dirname(parent.frame(2)$ofile) # get current dir

# manually enter the path (can automate it)
# also legend is rubbish, can use file.csv.legend to read the names of the categories
f <- file.path("/home/miroslav/Desktop/nmachinel.test/m05/TRAIN.D.og.op.nw.500_5.0_50.0_5/feat.csv")
# basic read of the csv file for visualiztions
train<-read.csv(f, header=F)
boxplot(V6~V10, data=train)
grid()

plot(train$V7, train$V8, asp=1)
