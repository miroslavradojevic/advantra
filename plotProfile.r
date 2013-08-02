rm(list = ls())

#cars <- c(1, 3, 6, 4, 9)
#plot(cars, type="l")

# CSV input path
filePath = "~/ImageJ/profile.csv"
#pdfPath  = "~/ImageJ/figure.pdf";

# read data
dat <- read.csv(file=filePath, head=TRUE, sep=",")

# Define colors
#plot_colors <- c(rgb(r=0.0,g=0.0,b=0.9), "red", "forestgreen")

# Start PDF device driver to save output to figure.pdf
#pdf(file=pdfPath, height=3.5, width=5)

# Trim off excess margin space (bottom, left, top, right)
par(mar=c(2.2, 2.0, 0.2, 0.2))

# Graph autos using a y axis that uses the full range of value
# in autos_data. Label axes with smaller font and use larger 
# line widths.
plot(dat$Angle, dat$Response, type="l", ylim=range(dat$Response), 
     axes=T, ann=F, xlab="Angle", ylab="Response", las=1, 
     cex.axis=0.8, lwd=2)

# Create box around plot
box()

# add grid
grid()