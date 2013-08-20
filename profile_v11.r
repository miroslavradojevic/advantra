rm(list = ls())

# input path
fileName="~/ImageJ/profile_v11.log"

conn=file(fileName, open="r")
linn=readLines(conn)

# PLOT PROFILE A
# line 1 profile angles in degrees
# line 2 profile values at ring A
# line 3 peak angles for A
# line 4 peak values for A

# export pdf
pdf(file='~/ImageJ/profileA.pdf', width = 20, height = 8) # in cm

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(4.2, 4.0, 0.5, 0.5))

# line 1
lineNr = 1; t <- unlist(strsplit(linn[lineNr], split=" ")); xName = t[1]
angA <- vector()
for (j in 2:length(t)){
  angA <- c(angA, as.numeric(t[j]))  
}

# line 2
lineNr = 2; t <- unlist(strsplit(linn[lineNr], split=" ")); yName = t[1]
profileA <- vector()
for (j in 2:length(t)){
  profileA <- c(profileA, as.numeric(t[j]))  
}

# plot this line
plot(angA, profileA, xlab=xName, ylab=yName, col="black", type="l", lwd=1, cex.axis=1.5) # , main="my_plot"
box()
grid()
lines(angA, profileA, xlab=xName, ylab=yName, col="black", type="o", lwd=1)
minProfileA = min(profileA)

# line 3,4 peaks A
lineNr = 3;
t <- unlist(strsplit(linn[lineNr], split=" "))
peakAnglesA <- vector()
for (j in 2:length(t)){
  peakAnglesA <- c(peakAnglesA, as.numeric(t[j]))  
}

lineNr = 4;
t <- unlist(strsplit(linn[lineNr], split=" "))
peakValuesA <- vector()
for (j in 2:length(t)){
  peakValuesA <- c(peakValuesA, as.numeric(t[j]))  
}

for (i in 1:length(peakValuesA)){
    segments(peakAnglesA[i], minProfileA, peakAnglesA[i], peakValuesA[i], col="red", lwd="4")  
}

dev.off()



close(conn)