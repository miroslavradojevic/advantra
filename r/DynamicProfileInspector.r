rm(list = ls())

# input path
fileName="~/DynamicProfileInspector.dat"


conn=file(fileName, open="r")
linn=readLines(conn)
# contains membership functions that fuzzify input variables and the decision
# results are fuzzy sets for each
# line 1 ang
# line 2 profile


# line 3 startX
# line 4 startY
# line 5 finishX
# line 6 finishY

# export pdf
pdf(file='~/profile.pdf', width = 8, height = 2) # in inches
# Trim off excess margin space (bottom, left, top, right)
par(mar=c(3.4, 3.0, 0.5, 0.5),mgp=c(2.4, 0.7, 0))

lineNr = 1;  # line 1
t <- unlist(strsplit(linn[lineNr], split=" "));
angA <- vector()
for (j in 2:length(t)){
  angA <- c(angA, as.numeric(t[j]))  
}

lineNr = 2;  # line 2
t <- unlist(strsplit(linn[lineNr], split=" ")); 
profile <- vector()
for (j in 2:length(t)){
  profile <- c(profile, as.numeric(t[j]))  
}

lineNr = 3;  # line 3
t <- unlist(strsplit(linn[lineNr], split=" ")); 
startX <- vector()
for (j in 2:length(t)){
  startX <- c(startX, as.numeric(t[j]))  
}

lineNr = 4;  # line 4
t <- unlist(strsplit(linn[lineNr], split=" ")); 
startY <- vector()
for (j in 2:length(t)){
  startY <- c(startY, as.numeric(t[j]))  
}


lineNr = 5;  # line 5
t <- unlist(strsplit(linn[lineNr], split=" ")); 
finishX <- vector()
for (j in 2:length(t)){
  finishX <- c(finishX, as.numeric(t[j]))  
}

lineNr = 6;  # line 6
t <- unlist(strsplit(linn[lineNr], split=" ")); 
finishY <- vector()
for (j in 2:length(t)){
  finishY <- c(finishY, as.numeric(t[j]))  
}

colors=c("black", "red", "green")

plot(angA, profile, 
     xlab=expression(alpha), 
     #axes=T,
     ann=T,
     ylab="",
     #ylim=c(50, 60),
     col=colors[1], 
     type="l", 
     lwd=3, 
     #las=1, ,
     #cex.axis=1.5,
     cex.lab=2.1,
     axes=F,
     frame.plot=T
)
axis(2, at=c(floor(min(profile)), ceiling(0.95*max(profile))), las=1, tck=-0.03, cex.axis=1.5)
axis(1, cex.axis=1.5,tck=-0.03)

points(finishX, finishY, 
       #xlab="", 
       #ylab="", 
       col=colors[2], 
       #type="o", 
       lwd=1,
       pch=4,
       cex=3.0
)

# points(startX, startY, 
#       xlab="", 
#       ylab="", 
#       col=colors[2], 
#       #type="o", 
#       lwd=2,
#       pch=20
# )

#legend("topright", c("PROFILE", "PEAKS"), fill= colors) ##, "START"
box()
grid()

dev.off()