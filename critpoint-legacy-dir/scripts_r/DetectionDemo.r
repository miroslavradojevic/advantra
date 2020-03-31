# Loading table from log file, space separated
# Create Grouped & Stacked Bar Plot

rm(list = ls())

# set working directory to current source directory
this.dir <- dirname(parent.frame(2)$ofile)
setwd(this.dir)

print(this.dir)

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
L = length(every1) # nr of folders

# add last row,sum per each column
tableDataCum = colSums (tableData, na.rm = FALSE, dims = 1)

tableData = rbind(tableData, tableDataCum)

# fill the values for output with zeros
plotData <- matrix(data=0, nrow=3,ncol=2*L)  # plotData for barplot

for (imageIdx in 1:(nrImages+1)) {

  print(imageIdx)

  for (i in 1:L) {

    a = tableData[imageIdx, every1[i]] # TP
    b = tableData[imageIdx, every2[i]] # FP
    c = tableData[imageIdx, every3[i]] # FN

    #print(c(a, b, c))

    # to normalize
    plotData[1,2*i-1] = a/(a+c)
    plotData[2,2*i-1] = c/(a+c)
    plotData[3,2*i  ] = b/(a+c)

    #plotData[1,2*i-1] = a
    #plotData[2,2*i-1] = c
    #plotData[3,2*i  ] = b

  }

  if (imageIdx<=nrImages) {
    filename = paste("ratesImg",imageIdx,".pdf", sep="")
  }
  else {
    filename = paste("alltogether.pdf", sep="")

    # in case of the last row (sum) add the values to plotDataTable so that each column corresponds to one folder (category)


  }

  pdf(file=filename, width = 5, height = 2) # in inches
  # Trim off excess margin space (bottom, left, top, right)
  # xpd=TRUE enable things to be drawn outside the plot region
  par(mar=c(1.9, 3.5, 0.5, 1.5),mgp=c(2.4, 0.7, 0), xpd=TRUE)

  colors=c(rgb(0,0,1,0.3),rgb(1,0,0,0.3),rgb(0,1,0,0.3))
  bardensity <- t(t(c(15,15,0)))
  barangle <- t(t(c(45,135,0)))
  spacing <-c(.75, .15)

  barplot(plotData,
          space=spacing,
          #xlab = "CONFIGURATION",
          ylab = "TP/FP/FN RATE",
          names.arg = legendData2,
          #legend.text = c("TP", "FP","FN"),
          col = colors,
          beside=FALSE

  )
  barplot(plotData,
          space=spacing,
          #xlab = "CONFIGURATION",
          ylab = "TP/FP/FN RATE",
          names.arg = legendData2,
          #legend.text = c("TP", "FP","FN"),
          col = "black",
          beside=FALSE,
          density = bardensity,
          angle = barangle,
          add=T
  )

  # Add legend to top right, outside plot region
  legend("topright", inset=c(-0.08,0), c("TP", "FN","FP"), fill= colors, bty = "n")
  dev.off()



}
# print last row for report table every 3rd goes to one column
print("table data")
print(tableData)
plotDataTable <- matrix(data=0, nrow=3,ncol=2*L) # plot in the table

for (i in 1:L) {

  a = tableData[nrImages+1, every1[i]] # TP
  b = tableData[nrImages+1, every2[i]] # FP
  c = tableData[nrImages+1, every3[i]] # FN

  plotDataTable[1,2*i-1] = a;
  plotDataTable[2,2*i-1] = b;
  plotDataTable[3,2*i-1] = c;

      plotDataTable[1,2*i] = a/(a+c);
      plotDataTable[2,2*i] = b/(a+c);
      plotDataTable[3,2*i] = c/(a+c);

}

print("table data TP (TPR);FP (FPR);FN in rows")
print(legendData2)
print(plotDataTable)

