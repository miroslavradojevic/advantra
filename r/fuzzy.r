rm(list = ls())

# input path
fileName="~/ImageJ/fuzzy.log"

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
pdf(file='~/ImageJ/h_LOW.pdf', width = 20, height = 10) # in cm

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(5.0, 6.5, 0.5, 0.5))

plot(theta, h_low, xlab=theta_name, ylab=h_low_name, 
     col="red", type="l", lwd=6, las=1, cex.axis=1.5, cex.lab=2.8, font = 1, family = "serif")
box(); grid()
dev.off()

#########  export pdf plot 3 vs. 1
pdf(file='~/ImageJ/h_MID.pdf', width = 20, height = 10) # in cm

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(5.0, 6.5, 0.5, 0.5))

plot(theta, h_mid, xlab=theta_name, ylab=h_mid_name, 
     col="blue", type="l", lwd=6, las=1, cex.axis=1.5, cex.lab=2.8, font = 1, family = "serif")
box(); grid()
dev.off()

#########  export pdf plot 4 vs. 1
pdf(file='~/ImageJ/h_HGH.pdf', width = 20, height = 10) # in cm

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(5.0, 6.5, 0.5, 0.5))

plot(theta, h_hgh, xlab=theta_name, ylab=h_hgh_name, 
     col="blue", type="l", lwd=6, las=1, cex.axis=1.5, cex.lab=2.8, font = 1, family = "serif")
box(); grid()
dev.off()

#########  export pdf plot 2,3,4 vs. 1
pdf(file='~/ImageJ/h.pdf', width = 20, height = 10) # in cm

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(5.0, 6.5, 0.5, 0.5))

plot(theta, h_hgh, xlab=theta_name, ylab=expression("h("*theta*")"), 
     col="red", type="l", lwd=6, las=1, cex.axis=1.5, cex.lab=2.8, font = 1, family = "serif")

lines(theta, h_mid, xlab=theta_name, ylab=h_mid_name, col="green", type="l", lwd=5)
lines(theta, h_low, xlab=theta_name, ylab=h_low_name, col="blue", type="l", lwd=4)
legend(-20,0.8, c("HIGH","MID","LOW"), lty=c(1,1,1), lwd=c(4,4,4),col=c("red","green","blue"), cex=2.5) 
# places a legend at the appropriate place, puts text in the legend, gives the legend appropriate symbols (lines)
# gives the legend lines the correct color and width
box(); grid()
dev.off()

#########  export pdf plot 6,7,8 vs. 5
pdf(file='~/ImageJ/q.pdf', width = 20, height = 10) # in cm

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(5.0, 6.5, 0.5, 0.5))

plot(x, q_yes, xlab=x_name, ylab=expression("q("*x*")"), 
     col="red", type="l", lwd=6, las=1, cex.axis=1.5, cex.lab=2.8, font = 1, family = "serif")
lines(x, q_maybe, xlab=x_name, ylab=q_maybe_name, col="green", type="l", lwd=5)
lines(x, q_no, xlab=x_name, ylab=q_no_name, col="blue", type="l", lwd=4)
legend(0,0.6, c("YES","MAY","NO"), lty=c(1,1,1), lwd=c(4,4,4),col=c("red","green","blue"), cex=2.5)
box(); grid()
dev.off()