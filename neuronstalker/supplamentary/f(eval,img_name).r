## clear 
rm(list = ls());

this.dir <- dirname(parent.frame(2)$ofile) # get current dir
source_dir <- file.path(this.dir, "summary") # directory to read input from

load(file=file.path(source_dir, "det.RData"))
# load(file=file.path(source_dir, "eval.RData"))
load(file=file.path(source_dir, "param.RData"))

df<-cbind(det, param) # , eval
remove(param, det ) # , eval

# no multiscale
#df<-df[df$D %in% c("4.0", "6.0", "8.0", "10.0"),]

# correct for the wrong tests
# df <- df[df$sthH %in% c("0.00", "5.00", "10.00"),]

################################################################################################
## take the one with min SD for each NAME (and bring along some more columns)
SD<-do.call(rbind, lapply(split(df, df[,"NAME"]), function(x) x[which.min(x[,"SD"]), ]) ) #c("NAME","eJUN", "D", "dsens", "nccL", "lhoodH", "sthH", "sthL")
cat("SD\n");print(SD)

pdf(file.path(source_dir, "f(SD,image).pdf"), width=4, height=3)
par(mar=c(1.2,3,0.5,0.5)+0.1, mgp=c(2.0,1,0)) 
mp<-barplot(SD$SD, ylab=expression('SD'), col=rgb(1,0,0,0.4),  xpd=F, las=2) #ylim=c(.0, 1),
text(mp, par("usr")[3] - 0.015, srt = 45, adj = 1, labels=SD$NAME, xpd = TRUE, font = 2, cex=0.8)
grid()
dev.off()

################################################################################################
## take the one with min SSD for each NAME (and bring along some more columns)
SSD<-do.call(rbind, lapply(split(df, df[,"NAME"]), function(x) x[which.min(x[,"SSD"]), ]) ) #c("NAME","eEND", "D", "dsens", "nccL", "lhoodH", "sthH", "sthL")
cat("SSD\n")
print(SSD)

pdf(file.path(source_dir, "f(SSD,image).pdf"), width=4, height=3)
par(mar=c(1.2,3,0.5,0.5)+0.1, mgp=c(2.0,1,0))
mp<-barplot(SSD$SSD, ylab=expression('SSD'), col=rgb(1,1,0,0.4), ylim=c(.0, 1), xpd=F, las=2)
text(mp, par("usr")[3] - 0.015, srt = 45, adj = 1, labels=SSD$NAME, xpd = TRUE, font = 2, cex=0.8)
grid()
dev.off()

################################################################################################
## take the one with min percSSD for each NAME (and bring along some more columns)
percSSD<-do.call(rbind, lapply(split(df, df[,"NAME"]), function(x) x[which.min(x[,"percSSD"]), ]) ) # c("NAME","eII", "eJUN", "eEND", "D", "dsens", "nccL", "lhoodH", "sthH", "sthL")
cat("percSSD\n")
print(percSSD)

pdf(file.path(source_dir, "f(percSSD,image).pdf"), width=8, height=3)
par(mar=c(1.2,3,0.5,0.5)+0.1, mgp=c(2.0,1,0))
mp<-barplot(percSSD$percSSD, ylab=expression('percSSD'), col=rgb(0,0,0,0.4), ylim=c(.0, 1), xpd=F, las=2)
text(mp, par("usr")[3] - 0.015, srt = 45, adj = 1, labels=percSSD$NAME, xpd = TRUE, font = 2, cex=0.8)
grid()
dev.off()

################################################################################################
## save data.frame out with [name,eI,eII,eJUN,eEND] to be used later on when concatenating data.frames for different snr-s
out<-cbind(
  data.frame(NAME=SD$NAME),
  data.frame(SD=SD[,"SD"]),
  data.frame(SSD=SSD[,"SSD"]), 
  data.frame(percSSD=percSSD[,"percSSD"])
)

save(out, file=file.path(source_dir, "out.RData"))

cat("TABLE (ISBI15) \n");
tabl <- rbind(as.character(SD$NAME), 
              format(SD$SD, digits=4), 
              format(SSD$SSD, digits=4), 
              format(percSSD$percSSD, digits=4)
              )
write.table(tabl, "", row.names=F, col.names=F, sep=" & ", eol=" \\\\ \n", quote=F)

stop("HEY!!!! OK");

out <- rbind(
  data.frame(t( SD[ ,"SD"  ] )),
  data.frame(t( SSD[,"SSD"] )),
  data.frame(t( percSSD[,"percSSD"] ))
)

names(out) <- levels(SD$NAME)
row.names(out) <- c("SD", "SSD", "percSSD")



stop("HEY!!!! OK");
# print table with the detection scores and parameters for different (4) evaluation categories
# write.table(  format(out, digits=2),
#               quote=F, row.names=T, col.names=NA,sep=" & ", eol=" \\\\ \n",
#               file=file.path(source_dir, "all.txt")
# );

cat("OUT: \n")
print(out)

# pdf(file.path(source_dir, "f(f+fJUN+fEND,image).pdf"), width=6, height=4)
# par(mar=c(3,2.5,1.5,1)+0.1, mgp=c(2, 1, 0))
# bp<-barplot(t(as.matrix(outII[,c("eJUN", "eII", "eEND")])), beside = T, col=c(rgb(1,0,0,0.4), rgb(0,0,0,0.6), rgb(1,1,0,0.4)), 
#             density=c(NA,10,NA), axisnames=F, ylim=c(.0, 1), xpd=F, 
#             width=c(0.2,1,0.2))
# text(colMeans(bp), par("usr")[3] - 0.010, srt = 45, adj = 1, labels=names(out), xpd = TRUE, font = 2, cex=0.6)
# legend("topright", c(expression('F'['JUN']), "F", expression('F'['END'])), bty="n", fill=c(rgb(1,0,0,0.4), rgb(0,0,0,0.6), rgb(1,1,0,0.4)), 
#        density=c(NA,10,NA), horiz=F, xpd=T, inset=c(-0.05,-0.1)) #, lwd=c(4,4,4) 
# grid()
# dev.off()

# pdf(file.path(source_dir, "f_boxplot.pdf"), width=6, height=4)
# par(mar=c(3,2.5,1.5,1)+0.1, mgp=c(2, 1, 0))
# boxplot(outII$eII, ylim=c(.5, 1)) # ylab=expression('F'['']), cex=0.5,
# grid()
# dev.off()
cat("MEDIAN \n")
cat("eI :  ", median(outI$eI),     "\n")
cat("eII:  ", median(outII$eII),   "\n")
cat("eJUN: ", median(outJUN$eJUN), "\n")
cat("eEND: ", median(outEND$eEND), "\n")


