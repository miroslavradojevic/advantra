rm(list = ls())

D     = 4         # neuron diameter
sigma = D/4       # profile 
td    = seq(-D/2, D/2,by=0.01)
w     = exp(-(td*td) / (2*sigma*sigma)) 

#pdf(file="~/filter-cs.pdf", width = 5, height = 4)

par(mar=c(4.2, 4.2, 1.0, 1.0))

plot(td, w, 
     type="l", 
     ylim=range(w), 
     #axes=T, 
     #ann=T, 
     xlab="d [pixels]", 
     #ylab="weight", 
     las=1, 
     cex.axis=1.8,
     cex.lab=2.5, 
     lwd=3
     )
box()
grid()
#dev.off()
