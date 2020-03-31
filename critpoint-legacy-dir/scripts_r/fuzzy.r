rm(list = ls())

# input path
fileName="~/fuzzy.dat"

conn=file(fileName, open="r")
linn=readLines(conn)

# contains membership functions that fuzzify input variables and the decision
# results are fuzzy sets
# line 1 theta
# line 2 h_low[theta] 
# line 3 h_mid[theta]
# line 4 h_hgh[theta]
# line 5 x ranging [0,1]
# line 6 Q_YES
# line 7 Q_MAYBE
# line 8 Q_NO

######### start reading ######### 
# load data from log file
# line 1 theta
lineNr = 1; t <- unlist(strsplit(linn[lineNr], split=" ")); 
theta_name = expression("" * theta)
theta <- vector()
for (j in 2:length(t)){
  theta <- c(theta, as.numeric(t[j]))  
}

# line 2 h_low[theta] 
lineNr = 2; t <- unlist(strsplit(linn[lineNr], split=" ")); 
h_low_name = expression("h LOW(" * theta * ")")
h_low <- vector()
for (j in 2:length(t)){
  h_low <- c(h_low, as.numeric(t[j]))  
}

# line 3 h_mid[theta] 
lineNr = 3; t <- unlist(strsplit(linn[lineNr], split=" "));
h_mid_name = expression("h MID(" * theta * ")")
h_mid <- vector()
for (j in 2:length(t)){
  h_mid <- c(h_mid, as.numeric(t[j]))  
}

# line 4 h_hgh[theta] 
lineNr = 4; t <- unlist(strsplit(linn[lineNr], split=" ")); 
h_hgh_name = expression("h HGH(" * theta * ")")
h_hgh <- vector()
for (j in 2:length(t)){
  h_hgh <- c(h_hgh, as.numeric(t[j]))  
}

# line 5 x 
lineNr = 5; t <- unlist(strsplit(linn[lineNr], split=" ")); 
x_name = expression("" * x * "")
x <- vector()
for (j in 2:length(t)){
  x <- c(x, as.numeric(t[j]))  
}

# line 6 q_yes[x] 
lineNr = 6; t <- unlist(strsplit(linn[lineNr], split=" ")); 
q_yes_name = expression("q YES(" * x * ")")
q_yes <- vector()
for (j in 2:length(t)){
  q_yes <- c(q_yes, as.numeric(t[j]))  
}

# line 7 q_maybe[x] 
lineNr = 7; t <- unlist(strsplit(linn[lineNr], split=" ")); 
q_maybe_name = expression("q MAYBE(" * x * ")")
q_maybe <- vector()
for (j in 2:length(t)){
  q_maybe <- c(q_maybe, as.numeric(t[j]))  
}

# line 8 q_no[x] 
lineNr = 8; t <- unlist(strsplit(linn[lineNr], split=" ")); 
q_no_name = expression("q MAYBE(" * x * ")")
q_no <- vector()
for (j in 2:length(t)){
  q_no <- c(q_no, as.numeric(t[j]))  
}

close(conn)
######### finished reading #########  

#########  export pdf plot 2 vs. 1
pdf(file='~/h_low.pdf', width = 10, height = 7)

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(5.0, 6.5, 0.5, 0.5))

plot(theta, h_low, 
     xlab=theta_name, 
     ylab=h_low_name, 
     col="red", 
     type="l", 
     lwd=6, 
     las=1, 
     #font = 1, 
     #family = "serif",
     cex.axis=1.8, 
     cex.lab=2.8 
)
box(); grid()
dev.off()

#########  export pdf plot 3 vs. 1
pdf(file='~/h_mid.pdf', width = 10, height = 7)

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(5.0, 6.5, 0.5, 0.5))

plot(theta, h_mid, 
     xlab=theta_name, 
     ylab=h_mid_name, 
     col="blue", 
     type="l", 
     lwd=6, 
     las=1, 
#      font = 1, 
#      family = "serif",  
     cex.axis=1.8, 
     cex.lab=2.8 
)
box(); grid()
dev.off()

#########  export pdf plot 4 vs. 1
pdf(file='~/h_hgh.pdf', width = 10, height = 7)

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(5.0, 6.5, 0.5, 0.5))

plot(theta, h_hgh, 
     xlab=theta_name, 
     ylab=h_hgh_name, 
     col="blue", 
     type="l", 
     lwd=6, 
     las=1, 
     #font = 1, 
     #family = "serif",
     cex.axis=1.8, 
     cex.lab=2.8 
)
box(); grid()
dev.off()


#################################
colors=c(rgb(1,0,0,1),rgb(0,1,0,1),rgb(0,0,1,1))

#########  export pdf plot 2,3,4 vs. 1
pdf(file='~/h.pdf', width = 5, height = 3)

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(3.3, 4.5, 0.5, 0.5),mgp=c(2.4, 0.7, 0))

plot(theta, h_hgh, 
     xlab=theta_name, 
     ylab=expression("h("*theta*")"), 
     col=colors[1], 
     type="l", 
     lwd=5, 
     las=1, 
     #font = 1, 
     #family = "serif",
     cex.axis=1.5, 
     cex.lab=2.1 ,
     lty=5,
     axes=F,
     frame.plot=T
)

axis(2, at=c(0, 0.5, 1.0),las=1,tck=-0.03,cex.axis=1.5)
axis(1,,cex.axis=1.5,tck=-0.03)
# points(theta, h_hgh, 
#        xlab="", 
#        ylab="", 
#        col="black", 
#        #type="o", 
#        lwd=2,
#        pch=20
# )

lines(theta, h_mid, 
      xlab=theta_name, 
      ylab=h_mid_name, 
      col=colors[2], 
      type="l", 
      lwd=5,
      lty=1
)
# points(theta, h_mid, 
#        xlab="", 
#        ylab="", 
#        col="black", 
#        #type="o", 
#        lwd=2,
#        pch=5
# )

lines(theta, h_low, 
      xlab=theta_name, 
      ylab=h_low_name,
      col=colors[3], 
      type="l", 
      lwd=5,
      lty=2
)
# points(theta, h_low, 
#        xlab="", 
#        ylab="", 
#        col="black", 
#        #type="o", 
#        lwd=2,
#        pch=12
# )

# plot on top again
lines(theta, h_hgh, 
     xlab=theta_name, 
     ylab=expression("h("*theta*")"), 
     col=colors[1], 
     type="l", 
     lwd=5, 
     las=1, 
     #font = 1, 
     #family = "serif",
     cex.axis=1.5, 
     cex.lab=2.1 ,
     lty=5
)

legend("left", 
       c("HIGH","MID","LOW"), 
       lty=c(1,1,1), 
       lwd=c(4,4,4),
       col=colors,
       bty = "n",
       cex=1.3
)  #, 
# places a legend at the appropriate place, puts text in the legend, gives the legend appropriate symbols (lines)
# gives the legend lines the correct color and width
box(); 
#grid()
dev.off()

#########  export pdf plot 6,7,8 vs. 5
pdf(file='~/q.pdf', width = 5, height = 3)

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(3.3, 4.5, 0.5, 0.5),mgp=c(2.4, 0.7, 0))

plot(x, q_yes, 
     xlab=x_name, 
     ylab=expression("q("*x*")"), 
     col=colors[1], 
     type="l", 
     lwd=5, 
     las=1, 
     #font = 1, 
     #family = "serif",
     cex.axis=1.5, 
     cex.lab=2.1,
     axes=F,
     frame.plot=T,
     lty=5
)
axis(2, at=c(0, 0.5, 1.0),las=1,tck=-0.03,cex.axis=1.5)
axis(1,,cex.axis=1.5,tck=-0.03)

lines(x, q_maybe, 
      xlab=x_name, 
      ylab=q_maybe_name, 
      col=colors[2], 
      type="l", 
      lwd=5,
      lty=1
)
lines(x, q_no, 
      xlab=x_name, 
      ylab=q_no_name, 
      col=colors[3], 
      type="l", 
      lwd=5,
      lty=2
)
legend("center", 
       c("YES","MAY","NO"), 
       lty=c(1,1,1), 
       lwd=c(4,4,4),
       col=colors,
       cex=1.3,
       bty = "n"
) # , 
#box(); 
#grid()
dev.off()