## clear 
rm(list = ls());

## functions 
eval.GENERAL    <- function(TP, FP, FN) {
  
  P <- ifelse(TP+FP>0, TP/(TP+FP), 0) 
  R <- ifelse(TP+FN>0, TP/(TP+FN), 0)
  
  eval_score <- ifelse(P+R>0, (2*P*R)/(P+R), 0)
  
  return(eval_score)
}

eval.I  <- function(TPj, FPj, FNj, TPe, FPe, FNe) {
  
  Pj <- ifelse(TPj+FPj>0, TPj/(TPj+FPj), 0)
  Rj <- ifelse(TPj+FNj>0, TPj/(TPj+FNj), 0)
  
  Pe <- ifelse(TPe+FPe>0, TPe/(TPe+FPe), 0)
  Re <- ifelse(TPe+FNe>0, TPe/(TPe+FNe), 0)
  
  Fj <- ifelse(Pj+Rj>0, (2*Pj*Rj)/(Pj+Rj), 0)
  Fe <- ifelse(Pe+Re>0, (2*Pe*Re)/(Pe+Re), 0)
  
  eval_1_score <- ifelse(Fj+Fe>0, (2*Fj*Fe)/(Fj+Fe), 0) 
  
  return(eval_1_score)
}

eval.II  <- function(TPj, FPj, FNj, TPe, FPe, FNe) {
  
  TP <- TPj + TPe
  FP <- FPj + FPe
  FN <- FNj + FNe
  
  P <- ifelse(TP+FP>0, TP/(TP+FP), 0)
  R <- ifelse(TP+FN>0, TP/(TP+FN), 0)
  
  eval_2_score <- ifelse(P+R>0, (2*P*R)/(P+R), 0)
  
  return(eval_2_score)
}

count.lines.csv <- function(file_path) {
  return(length(readLines(file_path)))
}

param.values <- function(file_path) {
  dir_names <- unlist(strsplit(file_path, "/"))
  # assume .csv was in the folder whose name carries the parameters (1 means one directory above)
  detection_dir_name <- dir_names[length(dir_names)-1] 
  detection_dir_name_parts <- unlist(strsplit(detection_dir_name, "_"))
  legend <- detection_dir_name_parts[1]
  legend_names <- unlist(strsplit(legend, "[.]"))
  tag <- legend_names[2:length(legend_names)]
  vals <- (detection_dir_name_parts[2:length(detection_dir_name_parts)])
  out <- data.frame()
  out <- rbind(out, vals)
  names(out) <- tag
  return(out)
}

this.dir <- dirname(parent.frame(2)$ofile) # get current dir
out_dir <- file.path(this.dir, "summary") # output directory for the script
dir.create(out_dir, showWarnings = FALSE)

cat("reading folders... ")

t1 <- proc.time()[1]

list = list.dirs(this.dir);  # will recursively list all the directories (with the detections inside)
if (length(list)<=0) stop("no folders with detections found");  

pamameter_records       <- list.files(list, pattern="eval.csv", full.names=T)

##
parameter_group_size    <- as.numeric(lapply(pamameter_records, count.lines.csv)) - 1 # first line is the legend, don't count it

##
params1 <- do.call("rbind", lapply(pamameter_records, param.values)) # extract parameters
param<-params1[rep(seq_len(nrow(params1)), parameter_group_size), ] # replicate with respect to the amount of lines in csv
remove(params1)

##
det<-do.call("rbind", lapply(pamameter_records, read.csv, header = TRUE, strip.white=TRUE))

# crop those that you don't need
drops <- c();#c("P_BIF","R_BIF","TP_BIF","FP_BIF","FN_BIF","P_CRS","R_CRS","TP_CRS","FP_CRS","FN_CRS")
det <- det[,!(names(det) %in% drops)]

t2 <- proc.time()[1] - t1
cat("DONE! elapsed ", t2, " seconds.\n")

# evaluation (det contains evaluation already)
# cat("evaluation... \n")
# eI <- eval.I(det$TP_JUN, det$FP_JUN, det$FN_JUN, det$TP_END, det$FP_END, det$FN_END)
# eII <- eval.II(det$TP_JUN, det$FP_JUN, det$FN_JUN, det$TP_END, det$FP_END, det$FN_END)
# eJUN <- eval.GENERAL(det$TP_JUN, det$FP_JUN, det$FN_JUN)
# eEND <- eval.GENERAL(det$TP_END, det$FP_END, det$FN_END)
# eval<-data.frame(eI=eI, eII=eII, eJUN=eJUN, eEND=eEND)

# save
cat("saving... \n")
save(det, file=file.path(out_dir, "det.RData"))
# save(eval, file=file.path(out_dir, "eval.RData"))
save(param, file=file.path(out_dir, "param.RData"))

cat("FINISHED")