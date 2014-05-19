## clear 
rm(list = ls())

ALLOC=100;
root_folder <- "/home/miroslav/Desktop/varyp.snr.dmin.nrimg._3.0_3.0_1/";   # choose the folder to summarize scores
list = list.dirs(root_folder);  # will recursively list everything
print(list);

leg <- c("scale", "Dlist");
det <- data.frame(
  Dmax=numeric(ALLOC),
  SNR=numeric(ALLOC),
  scale_bch=numeric(ALLOC), 
  p1=numeric(ALLOC),
  p2=numeric(ALLOC),
  p3=numeric(ALLOC),
  s=numeric(ALLOC),
  Dlist=character(ALLOC),
  M=numeric(ALLOC),
  L=numeric(ALLOC),
  nccH=numeric(ALLOC),
  nccL=numeric(ALLOC),
  lhoodH=numeric(ALLOC),
  lhoodL=numeric(ALLOC),
  outSig=numeric(ALLOC),
  TP_BIF=numeric(ALLOC),
  FP_BIF=numeric(ALLOC),
  FN_BIF=numeric(ALLOC),
  P_BIF=numeric(ALLOC),
  R_BIF=numeric(ALLOC),
  TP_END=numeric(ALLOC),
  FP_END=numeric(ALLOC),
  FN_END=numeric(ALLOC),
  P_END=numeric(ALLOC),
  R_END=numeric(ALLOC),
  stringsAsFactors = FALSE
  )

counter <- 0;

if (length(list)>0) {

  for (i in 1:length(list)) {
    
    if (grepl("synjun", list[i]) && grepl("det", list[i])) {
      
      #print(list[i]);
      
      elems <- strsplit(list[i], "/");
      
      # take out those elements that contain keywords
      found_synjun <- F;
      found_det    <- F;
      found_log    <- F;
      
      for (j in 1:length(elems[[1]])) {
        
        if (grepl("synjun", elems[[1]][j])) { # folder with generated junctions
          
          gen_params <- strsplit(elems[[1]][j], "_");
          
          # extract parameters
          Dmax      <- as.numeric(gen_params[[1]][2]); 
          SNR       <- as.numeric(gen_params[[1]][3]);          
          scale_bch <- as.numeric(gen_params[[1]][4]);
          p1        <- as.numeric(gen_params[[1]][5]);
          p2        <- as.numeric(gen_params[[1]][6]);
          p3        <- as.numeric(gen_params[[1]][7]);
          
          found_synjun <- T; 
          
        }
        
        
        
        if (grepl("det", elems[[1]][j])) { # folder with generated detections
          
          det_params <- strsplit(elems[[1]][j], "_");
          
          # extract parameters
          s       <- as.numeric(det_params[[1]][2]);
          Dlist   <-            det_params[[1]][3];
          M       <- as.numeric(det_params[[1]][4]);
          L       <- as.numeric(det_params[[1]][5]);
          nccH    <- as.numeric(det_params[[1]][6]);
          nccL    <- as.numeric(det_params[[1]][7]);
          lhoodH  <- as.numeric(det_params[[1]][8]);
          lhoodL  <- as.numeric(det_params[[1]][9]);
          outSig  <- as.numeric(det_params[[1]][10]);
          
          found_det    <- T;
          
        }
        
      }
      
      list_detections <- list.files(list[i], pattern="det.csv", full.names=T); # list[i] is current folder
      
      if (length(list_detections)==1) {
        
        res <- read.csv(list_detections[1], header=T);
        
        tp_bif <- sum(res$TP_BIF);
        fp_bif <- sum(res$FP_BIF);
        fn_bif <- sum(res$FN_BIF);
        
        p_bif <- tp_bif / (tp_bif+fp_bif);
        r_bif <- tp_bif / (tp_bif+fn_bif);
        
        tp_end <- sum(res$TP_END);
        fp_end <- sum(res$FP_END);
        fn_end <- sum(res$FN_END);
        
        p_end <- tp_end / (tp_end+fp_end);
        r_end <- tp_end / (tp_end+fn_end);
        
        found_log <- T;
        
      }
      
      if (found_synjun && found_det && found_log) {
      
        counter <- counter + 1;
        
        det$Dmax[counter] <- Dmax;
        det$SNR[counter] <- SNR;
        det$scale_bch[counter] <- scale_bch;
        det$p1[counter] <- p1;
        det$p2[counter] <- p2;
        det$p3[counter] <- p3;
        
        det$s[counter] <- s;
        det$Dlist[counter] <- Dlist;
        det$M[counter] <- M;
        det$L[counter] <- L;
        det$nccH[counter] <- nccH;
        det$nccL[counter] <- nccL;
        det$lhoodH[counter] <- lhoodH;
        det$lhoodL[counter] <- lhoodL;
        det$outSig[counter] <- outSig;
        
        det$TP_BIF[counter] <- tp_bif;
        det$FP_BIF[counter] <- fp_bif;
        det$FN_BIF[counter] <- fn_bif;
        det$P_BIF[counter]  <- p_bif;
        det$R_BIF[counter]  <- r_bif;
        det$TP_END[counter] <- tp_end;
        det$FP_END[counter] <- fp_end;
        det$FN_END[counter] <- fn_end;
        det$P_END[counter]  <- p_end;
        det$R_END[counter]  <- r_end;
        
      }
      
    }
    
  }
  
  det <- det[1:counter,];
  print(det);
  save(det, file=paste(root_folder,"det.Rda"));
  
} else {
  print("nothing");
}