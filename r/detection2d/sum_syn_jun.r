root_folder <- "~/vary_p"; # choose the folder to summarize scores
list = list.dirs(root_folder); # will recursively list everything

leg <- c("scale", "Dlist");
det <- data.frame(s=numeric(10), Dlist=character(10), stringsAsFactors = FALSE)

counter <- 0;

for (i in 1:length(list)) {
  
  if (grepl("syn_jun", list[i]) && grepl("det", list[i])) {
    
    elems <- strsplit(list[i], "/");
     
    # take out those elements that contain keywords
    found1 <- false;
    for (j in 1:length(elems[[1]])) {
       
      if (grepl("syn_jun", elems[[1]][j])) {
        gen_params <- strsplit(elems[[1]][j], "_");
        found1 <- 
      }
       
      if (grepl("det", elems[[1]][j])) {
        det_params <- strsplit(elems[[1]][j], "_");
        s<-as.numeric(det_params[[1]][2]);
        Dlist<-det_params[[1]][3];
      }

    }
    
    #counter <- counter + 1;
    #print(counter);
    #print(paste(s, Dlist, sep="|"));
    det$s[counter] <- s;
    det$Dlist[counter] <- Dlist;
  
  }
  
}

det <- det[1:counter,];
str(det);     
  #if (cnt==2) {
  #print(list[i]);
  #}
  
#   data_dir <- paste(root_folder,list_0[i], sep="/");
#   print(data_dir);
#   
#   list_1 <- list.dirs(data_dir);
#   
#   for (j in 1:length(list_1)) {
#     
#     det_dir <- paste(root_folder,list_0[i],list_1[j], sep="/");
#     print(det_dir);
#     
#   }