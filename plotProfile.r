rm(list = ls())

#cars <- c(1, 3, 6, 4, 9)
#plot(cars, type="l")

# CSV input path
filePath = "~/ImageJ/profile.csv"
msPath  = "~/ImageJ/ms.csv";
csPath  = "~/ImageJ/cs.csv";

# read data
dat <- read.csv(file=filePath, head=TRUE, sep=",")
ms <- read.csv(file=msPath, head=TRUE, sep=",")
cs <- read.csv(file=csPath, head=TRUE, sep=",")
# Define colors
#plot_colors <- c(rgb(r=0.0,g=0.0,b=0.9), "red", "forestgreen")

# Start PDF device driver to save output to figure.pdf
#pdf(file=pdfPath, height=3.5, width=5)

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(4.2, 4.0, 0.2, 0.2))

# Graph autos using a y axis that uses the full range of value
# in autos_data. Label axes with smaller font and use larger 
# line widths.
plot(dat$Angle, dat$Response, type="l", ylim=range(dat$Response), 
     axes=T, ann=T, xlab="orientation[deg]", ylab="profile response", las=1, 
     cex.axis=0.8, lwd=2)

lines(ms$Angle, ms$Response, type="p", col="green", pch=18)

#segments(cs$Angle, cs$Lower, cs$Angle, cs$Higher, col="red", lwd="5") 

# Create box around plot
box()

# add grid
grid()