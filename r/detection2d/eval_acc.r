# expects detection data frame to be loaded
# stops if not loaded
# script will check how F measures depend on input parameters such as ncc_low, likelihood_high

k_height <- 0.05;

if (length(ls())==0) stop("No variables");

exists <- F;
for (i in 1:length(ls())) 
  if (ls()[i]=="det") exists <- T;

if (!exists) stop("No dataframe det");

# extract diameters used for the detection as numnerical
dd <- as.numeric(det$Dlist);
k <- dd / det$Dmax;
# take set of indexes where k has particular level
# using hierarchical clustering, define the height of each cluster to include
idxs <- cutree(hclust(dist(k)), h=k_height);

for (i in 1:max(idxs)) { # for each cluster (detection D)
  
  idxs1 <- which(idxs %in% i);
  #print(idxs1);
  values_det[idxs1,]$nccL 
  values_det[idxs1,]$lhoodH 
  F1[idxs1]
  # recompose the matrix wit so that it 
  
  
}

beta <- 1;
F1_0_BIF <- (1 + beta^2) * ((det$P_BIF * det$R_BIF) / (beta^2 * det$P_BIF + det$R_BIF));
beta <- 0.5;
F0_5_BIF <- (1 + beta^2) * ((det$P_BIF * det$R_BIF) / (beta^2 * det$P_BIF + det$R_BIF));
beta <- 2;
F2_0_BIF <- (1 + beta^2) * ((det$P_BIF * det$R_BIF) / (beta^2 * det$P_BIF + det$R_BIF));
G_BIF <- sqrt(det$P_BIF * det$R_BIF);

beta <- 1;
F1_0_END <- (1 + beta^2) * ((det$P_END * det$R_END) / (beta^2 * det$P_END + det$R_END));
beta <- 0.5;
F0_5_END <- (1 + beta^2) * ((det$P_END * det$R_END) / (beta^2 * det$P_END + det$R_END));
beta <- 2;
F2_0_END <- (1 + beta^2) * ((det$P_END * det$R_END) / (beta^2 * det$P_END + det$R_END));
G_END <- sqrt(det$P_END * det$R_END);





#pdf("ncc_low.pdf") ;
#plot(det$nccL, F, main="NCC dependence", xlab="ncc_low", ylab="F score");
#dev.off();
#par(mfrow=c(1,1));
# scatter2D(det$nccL, 
#           det$lhoodH, 
#           #det$Dmax,
#           colvar=array(F,c(length(F),1)),
#           col=gg.col(length(F)),
#           #theta=150, 
#           #phi=10, 
#           type="p", 
#           pch=20); 
# , asp=1


#plot(det$lhoodH, F, main="LHOOD dependence", xlab="ncc_low", ylab="F score");
#plot(det$lhoodH, F, main="k dependence", xlab="ncc_low", ylab="F score");
