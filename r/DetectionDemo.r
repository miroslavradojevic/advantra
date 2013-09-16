# Loading table from log file, space separated
# Create Grouped & Stacked Bar Plot

rm(list = ls())

# set working directory to current source directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)
#X11()
######## data #########
fileNameData = "data.dat"
tableData=read.table(fileNameData)
######## legend #########
fileNameLegend = "legend.dat"
legendData=read.table(fileNameLegend)

legendData2 <- c() #(data="",nrow=1,ncol=2*length(legendData))
for (name in names(legendData)) {
    legendData2 <- c(legendData2, as.character(as.character(legendData[[name]])), "")
}

nrElements  = dim(tableData)[2]
nrImages    = dim(tableData)[1]

# each block one SNR, with 7 categories for that SNR

# reorganize data matrix
loopAll <- 1:nrElements
every1 <- loopAll[seq(1, length(loopAll), 3)]
every2 <- loopAll[1+seq(1, length(loopAll), 3)]
every3 <- loopAll[2+seq(1, length(loopAll), 3)]
L = length(every1)

# add last row,sum per each column
tableDataCum = colSums (tableData, na.rm = FALSE, dims = 1)

tableData = rbind(tableData, tableDataCum)

# fill the values for output with zeros
plotData <- matrix(data=0,nrow=3,ncol=2*L)

for (imageIdx in 1:(nrImages+1)) {
  print(imageIdx)

  for (i in 1:L) {
    a = tableData[imageIdx, every1[i]]
    b = tableData[imageIdx, every2[i]]
    c = tableData[imageIdx, every3[i]]

    #print(c(a, b, c))

    # to normalize
    plotData[1,2*i-1] = a/(a+c)
    plotData[2,2*i-1] = c/(a+c)
    plotData[3,2*i  ] = b/(a+c)

    #plotData[1,2*i-1] = a
    #plotData[2,2*i-1] = b
    #plotData[3,2*i  ] = c

  }

  if (imageIdx<=nrImages) {
    filename = paste("ratesImg",imageIdx,".pdf", sep="")
  }
  else {
    filename = paste("alltogether.pdf", sep="")
  }

  pdf(file=filename, width = 20, height = 10) # in cm
  colors=c(rgb(0,0,1,0.3),rgb(1,0,0,0.3),rgb(0,1,0,0.3))
  bardensity <- t(t(c(15,15,0)))
  barangle <- t(t(c(45,135,0)))
  barplot(plotData,
          space=c(.75,.25),
          #xlab = "JUNCTION CONFIGURATION",
          ylab = "TP/FP/FN RATE",
          names.arg = legendData2,
          #legend.text = c("TP", "FP","FN"),
          col = colors,
          beside=FALSE
  )
  barplot(plotData,
          space=c(.75,.25),
          #xlab = "JUNCTION CONFIGURATION",
          ylab = "TP/FP/FN RATE",
          names.arg = legendData2,
          #legend.text = c("TP", "FP","FN"),
          col = "black",
          beside=FALSE,
	  density = bardensity,
	  angle = barangle,
          add=T
  )
  legend("topright", c("TP", "FN","FP"), fill= colors)
  dev.off()

}
