# color
# red
> colorRampPalette(brewer.pal(8, "YlOrRd"))(15)
 [1] "#FFFFCC" "#FFF6B6" ["#FFEDA0"] "#FEE38B" ["#FED976"] "#FEC561" "#FEB24C" ["#FD9F44"] "#FD8D3C" ["#FC6D33"] "#FC4E2A" "#EF3322" ["#E31A1C"] "#CA0D20" ["#B10026"]
# green
> colorRampPalette(brewer.pal(8, "Greens"))(15)
 [1] "#F7FCF5" "#EEF8EA" "#E5F5E0" "#D6EFD0" "#C7E9C0" "#B4E1AD" "#A1D99B" "#8ACE88" ["#74C476"] "#5AB769" "#41AB5D" "#319B50" "#238B45" ["#11723B"] "#005A32"



### QTL all BSA6

setwd("D:/Dropbox (TX Biomed)/Emily/4.P01/6.BSA6.Mal31xKH004/3.analysis/BSA")

library(magrittr)
library(dplyr)
library("QTLseqr")
library("ggpubr")

format_genomic <- function(...) {
      function(x) {
            limits <- c(1e0,   1e3, 1e6)
            #prefix <- c("","Kb","Mb")
            # Vector with array indices according to position in intervals
            i <- findInterval(abs(x), limits)
            # Set prefix to " " for very small values < 1e-24
            i <- ifelse(i==0, which(limits == 1e0), i)
            paste(format(round(x/limits[i], 1),
                         trim=TRUE, scientific=FALSE, ...)
                #  ,prefix[i]
            )
      }
}

BSA <- read.table("BSA6.1.multi.MAL31xKH004.SNP.filter.table", sep = '\t',header = TRUE)
> dim(BSA)
[1] 12803   460

AD <- BSA[,seq(5,460,4)]
DP <- BSA[,seq(6,460,4)]
write.csv(DP, file = "BSA6.1.ALL.DP.csv",row.names=FALSE)

ref.DP <- function(X){as.numeric(strsplit(as.character(X),",")[[1]])[1]}
dim(AD)
[1] 12803   114

refFre.AD <- matrix(ncol=114,nrow=12803)
for(j in 1:114){
	refFre.AD[,j]<-sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j])
	}
	
refFre.AD[DP<20]<-NA
colnames(refFre.AD)<-colnames(AD)	
refFre.AD <- cbind(BSA[,1:4],refFre.AD)
colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))


write.csv(refFre.AD, file = "BSA6.1.ALL.refFre.AD.csv",row.names=FALSE)

######################################################
######### remove outliers ############################
######################################################
outliersMAD <- function(data, MADCutOff = 2, replace = NA, values = FALSE, bConstant = 1.4826, digits = 2) {
    absMADAway <- abs(   (data - median(data, na.rm = T))  /  mad(data, constant = bConstant, na.rm = T)  )
   
    data[absMADAway > MADCutOff] <- replace
    
    if (values == TRUE) { 
        return(round(absMADAway, digits)) #if values == TRUE, return number of mads for each value
    } else {
        return(round(data, digits)) #otherwise, return values with outliers replaced
    }
}

outlierByMAD <- function (x, k){
	n <- length(x)
     y <- x
     
 for (i in (k + 1):(n - k)) {
 	
    data <- x[(i - k):(i + k)]
    y[i] <- outliersMAD(data)[k+1]}
    
   return(y)
                               }

C37.M11.D10.filter <- outlierByMAD(refFre.AD$C37.M11.D10,50)

refFre.AD.BSA6.LC <- cbind(refFre.AD[,1:4], C37.M11.D10.filter,C37.M11.D21.filter,C37.M11.D3.filter,C37.M21.D10.filter,C37.M21.D21.filter,C37.M21.D3.filter,C39.M11.D10.filter,C39.M11.D21.filter,C39.M11.D3.filter,C39.M21.D10.filter,C39.M21.D21.filter,C39.M21.D3.filter,D.M11.T10.filter,D.M11.T2.filter,D.M11.T5.filter,D.M11.T8.filter,D.M21.T10.filter,D.M21.T2.filter,D.M21.T5.filter,D.M21.T8.filter,FG.BC.0196.filter,FG.BC.0197.filter,FG.BC.0198.filter,FG.BC.0199.filter,FG.BC.0200.filter,FG.BC.0201.filter,FG.BC.0202.filter,FG.BC.0203.filter,FG.BC.0204.filter,FG.BC.0207.filter,FG.BC.0210.filter,FG.BC.0213.filter,FG.BC.0216.filter,FG.BC.0219.filter,FG.BC.0222.filter,FG.BC.0225.filter,FG.BC.0228.filter,FG.BC.0231.filter,FG.BC.0234.filter,FG.BC.0237.filter,FG.BC.0240.filter,FG.BC.0243.filter,FG.BC.0246.filter,FG.BC.0249.filter,FG.BC.0252.filter,FG.BC.0255.filter,FG.BC.0258.filter,FG.BC.0261.filter,FG.BC.0264.filter,FG.BC.0267.filter,FG.BC.0270.filter,FG.BC.0273.filter,FG.BC.0276.filter,FG.BC.0279.filter,FG.BC.0282.filter,FG.BC.0285.filter,FG.BC.0288.filter,FG.BC.0300.filter,FG.BC.0303.filter,FG.BC.0306.filter,FG.BC.0318.filter,FG.BC.0321.filter,FG.BC.0324.filter,FG.BC.0327.filter,FG.BC.0330.filter,FG.BC.0333.filter,FG.BC.0336.filter,FG.BC.0339.filter,FG.BC.0342.filter,FG.BC.0354.filter,FG.BC.0357.filter,FG.BC.0360.filter,FG.BSA.0007.filter,FG.BSA.0008.filter,FG.BSA.0015.filter,FG.BSA.0016.filter,FG.BSA.0017.filter,FG.BSA.0018.filter,FG.BSA.0043.filter,FG.BSA.0049.filter,FG.BSA.0067.filter,FG.BSA.0075.filter,FG.BSA.0115.filter,FG.BSA.0123.filter,FG.BSA.0163.filter,FG.BSA.0171.filter,FG.BSA.0211.filter,FG.BSA.0219.filter,FG.BSA.0259.filter,FG.BSA.0267.filter,M.M11.T10.filter,M.M11.T2.filter,M.M11.T5.filter,M.M11.T8.filter,M.M21.T10.filter,M.M21.T2.filter,M.M21.T5.filter,M.M21.T8.filter,N.M1.filter,N.M1.T0.2.filter,N.M11.T10.filter,N.M11.T2.filter,N.M11.T5.filter,N.M11.T8.filter,N.M2.filter,N.M2.T0.2.filter,N.M21.T10.filter,N.M21.T2.filter,N.M21.T5.filter,N.M21.T8.filter,N.M3.filter,N.M3.T0.2.filter,N.M4.filter,N.M4.T0.2.filter)

> dim(refFre.AD.BSA6.LC)
[1] 12803   118


write.csv(refFre.AD.BSA6.LC, file = "BSA6.ALL.refFre.AD.filter.csv",row.names=FALSE)


refFre.AD.BSA6.LC.save <- refFre.AD.BSA6.LC

################################################################  
###### compare whole genome allele frequency  ##################  
################################################################  
library(reshape2)
library(ggplot2)
library(ggridges)

refFre.LC <- setNames(melt(refFre.AD.BSA6.LC[,5:118]), c('BSAs', 'Allele frequency of Mal31'))

library(doBy)
colnames(refFre.LC)[2]<- c("AlleleFrequency")
refFre.LC.filter <- refFre.LC[rowSums(is.na(refFre.LC)) == 0,]
sum <- summaryBy(AlleleFrequency ~ BSAs, data = refFre.LC.filter, FUN = list(mean, median))

write.csv(sum, file = "Mal31-AF.sum_coreGenome.ALL.csv",row.names=FALSE)



################################################################
######################## Allele frequency ######################
################################################################

#######tricubesmooth#######		
tricubeStat<-function(POS,Stat,windowSize=0.5e5,...)
{
if(windowSize<=0)
stop("Apositivesmoothingwindowisrequired")
stats::predict(locfit::locfit(Stat~locfit::lp(POS,h=windowSize,deg=0),...),POS)
}

#### too low coverage to analysis
dplyr::mutate(D.M11.T5.tricube=tricubeStat(POS=POS,Stat=D.M11.T5.filter,windowSize,...))
dplyr::mutate(FG.BC.0237.tricube=tricubeStat(POS=POS,Stat=FG.BC.0237.filter,windowSize,...))
dplyr::mutate(FG.BC.0330.tricube=tricubeStat(POS=POS,Stat=FG.BC.0330.filter,windowSize,...))
dplyr::mutate(FG.BSA.0008.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0008.filter,windowSize,...))
dplyr::mutate(FG.BSA.0015.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0015.filter,windowSize,...))
dplyr::mutate(FG.BSA.0018.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0018.filter,windowSize,...))
dplyr::mutate(M.M11.T5.tricube=tricubeStat(POS=POS,Stat=M.M11.T5.filter,windowSize,...))
dplyr::mutate(M.M21.T2.tricube=tricubeStat(POS=POS,Stat=M.M21.T2.filter,windowSize,...))
dplyr::mutate(M.M21.T5.tricube=tricubeStat(POS=POS,Stat=M.M21.T5.filter,windowSize,...))
dplyr::mutate(N.M11.T2.tricube=tricubeStat(POS=POS,Stat=N.M11.T2.filter,windowSize,...))


refFre.AD.BSA6.LC <- refFre.AD.BSA6.LC.save

SAnalysis.1<-function(SNPset,windowSize=1e5,...)
{SNPset<-SNPset%>%
dplyr::group_by(CHROM)%>%
dplyr::mutate(C37.M11.D10.tricube=tricubeStat(POS=POS,Stat=C37.M11.D10.filter,windowSize,...))%>%
dplyr::mutate(C37.M11.D3.tricube=tricubeStat(POS=POS,Stat=C37.M11.D3.filter,windowSize,...))%>%
dplyr::mutate(C37.M11.D21.tricube=tricubeStat(POS=POS,Stat=C37.M11.D21.filter,windowSize,...))%>%
dplyr::mutate(C37.M21.D10.tricube=tricubeStat(POS=POS,Stat=C37.M21.D10.filter,windowSize,...))%>%
dplyr::mutate(C37.M21.D21.tricube=tricubeStat(POS=POS,Stat=C37.M21.D21.filter,windowSize,...))%>%
dplyr::mutate(C37.M21.D3.tricube=tricubeStat(POS=POS,Stat=C37.M21.D3.filter,windowSize,...))%>%
dplyr::mutate(C39.M11.D10.tricube=tricubeStat(POS=POS,Stat=C39.M11.D10.filter,windowSize,...))%>%
dplyr::mutate(C39.M11.D21.tricube=tricubeStat(POS=POS,Stat=C39.M11.D21.filter,windowSize,...))%>%
dplyr::mutate(C39.M11.D3.tricube=tricubeStat(POS=POS,Stat=C39.M11.D3.filter,windowSize,...))%>%
dplyr::mutate(C39.M21.D10.tricube=tricubeStat(POS=POS,Stat=C39.M21.D10.filter,windowSize,...))%>%
dplyr::mutate(C39.M21.D21.tricube=tricubeStat(POS=POS,Stat=C39.M21.D21.filter,windowSize,...))%>%
dplyr::mutate(C39.M21.D3.tricube=tricubeStat(POS=POS,Stat=C39.M21.D3.filter,windowSize,...))%>%
dplyr::mutate(D.M11.T10.tricube=tricubeStat(POS=POS,Stat=D.M11.T10.filter,windowSize,...))%>%
dplyr::mutate(D.M11.T2.tricube=tricubeStat(POS=POS,Stat=D.M11.T2.filter,windowSize,...))%>%
dplyr::mutate(D.M11.T8.tricube=tricubeStat(POS=POS,Stat=D.M11.T8.filter,windowSize,...))%>%
dplyr::mutate(D.M21.T10.tricube=tricubeStat(POS=POS,Stat=D.M21.T10.filter,windowSize,...))%>%
dplyr::mutate(D.M21.T2.tricube=tricubeStat(POS=POS,Stat=D.M21.T2.filter,windowSize,...))%>%
dplyr::mutate(D.M21.T5.tricube=tricubeStat(POS=POS,Stat=D.M21.T5.filter,windowSize,...))%>%
dplyr::mutate(D.M21.T8.tricube=tricubeStat(POS=POS,Stat=D.M21.T8.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0196.tricube=tricubeStat(POS=POS,Stat=FG.BC.0196.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0197.tricube=tricubeStat(POS=POS,Stat=FG.BC.0197.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0198.tricube=tricubeStat(POS=POS,Stat=FG.BC.0198.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0199.tricube=tricubeStat(POS=POS,Stat=FG.BC.0199.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0200.tricube=tricubeStat(POS=POS,Stat=FG.BC.0200.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0201.tricube=tricubeStat(POS=POS,Stat=FG.BC.0201.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0202.tricube=tricubeStat(POS=POS,Stat=FG.BC.0202.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0203.tricube=tricubeStat(POS=POS,Stat=FG.BC.0203.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0204.tricube=tricubeStat(POS=POS,Stat=FG.BC.0204.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0207.tricube=tricubeStat(POS=POS,Stat=FG.BC.0207.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0210.tricube=tricubeStat(POS=POS,Stat=FG.BC.0210.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0213.tricube=tricubeStat(POS=POS,Stat=FG.BC.0213.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0216.tricube=tricubeStat(POS=POS,Stat=FG.BC.0216.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0219.tricube=tricubeStat(POS=POS,Stat=FG.BC.0219.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0222.tricube=tricubeStat(POS=POS,Stat=FG.BC.0222.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0225.tricube=tricubeStat(POS=POS,Stat=FG.BC.0225.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0228.tricube=tricubeStat(POS=POS,Stat=FG.BC.0228.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0231.tricube=tricubeStat(POS=POS,Stat=FG.BC.0231.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0234.tricube=tricubeStat(POS=POS,Stat=FG.BC.0234.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0240.tricube=tricubeStat(POS=POS,Stat=FG.BC.0240.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0243.tricube=tricubeStat(POS=POS,Stat=FG.BC.0243.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0246.tricube=tricubeStat(POS=POS,Stat=FG.BC.0246.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0249.tricube=tricubeStat(POS=POS,Stat=FG.BC.0249.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0252.tricube=tricubeStat(POS=POS,Stat=FG.BC.0252.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0255.tricube=tricubeStat(POS=POS,Stat=FG.BC.0255.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0258.tricube=tricubeStat(POS=POS,Stat=FG.BC.0258.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0261.tricube=tricubeStat(POS=POS,Stat=FG.BC.0261.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0264.tricube=tricubeStat(POS=POS,Stat=FG.BC.0264.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0267.tricube=tricubeStat(POS=POS,Stat=FG.BC.0267.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0270.tricube=tricubeStat(POS=POS,Stat=FG.BC.0270.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0273.tricube=tricubeStat(POS=POS,Stat=FG.BC.0273.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0276.tricube=tricubeStat(POS=POS,Stat=FG.BC.0276.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0279.tricube=tricubeStat(POS=POS,Stat=FG.BC.0279.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0282.tricube=tricubeStat(POS=POS,Stat=FG.BC.0282.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0285.tricube=tricubeStat(POS=POS,Stat=FG.BC.0285.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0288.tricube=tricubeStat(POS=POS,Stat=FG.BC.0288.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0300.tricube=tricubeStat(POS=POS,Stat=FG.BC.0300.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0303.tricube=tricubeStat(POS=POS,Stat=FG.BC.0303.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0306.tricube=tricubeStat(POS=POS,Stat=FG.BC.0306.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0318.tricube=tricubeStat(POS=POS,Stat=FG.BC.0318.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0321.tricube=tricubeStat(POS=POS,Stat=FG.BC.0321.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0324.tricube=tricubeStat(POS=POS,Stat=FG.BC.0324.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0327.tricube=tricubeStat(POS=POS,Stat=FG.BC.0327.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0333.tricube=tricubeStat(POS=POS,Stat=FG.BC.0333.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0336.tricube=tricubeStat(POS=POS,Stat=FG.BC.0336.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0339.tricube=tricubeStat(POS=POS,Stat=FG.BC.0339.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0342.tricube=tricubeStat(POS=POS,Stat=FG.BC.0342.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0354.tricube=tricubeStat(POS=POS,Stat=FG.BC.0354.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0357.tricube=tricubeStat(POS=POS,Stat=FG.BC.0357.filter,windowSize,...))%>%
dplyr::mutate(FG.BC.0360.tricube=tricubeStat(POS=POS,Stat=FG.BC.0360.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0007.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0007.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0016.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0016.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0017.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0017.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0043.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0043.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0049.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0049.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0067.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0067.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0075.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0075.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0115.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0115.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0123.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0123.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0163.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0163.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0171.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0171.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0211.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0211.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0219.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0219.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0259.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0259.filter,windowSize,...))%>%
dplyr::mutate(FG.BSA.0267.tricube=tricubeStat(POS=POS,Stat=FG.BSA.0267.filter,windowSize,...))%>%
dplyr::mutate(M.M11.T10.tricube=tricubeStat(POS=POS,Stat=M.M11.T10.filter,windowSize,...))%>%
dplyr::mutate(M.M11.T8.tricube=tricubeStat(POS=POS,Stat=M.M11.T8.filter,windowSize,...))%>%
dplyr::mutate(M.M21.T10.tricube=tricubeStat(POS=POS,Stat=M.M21.T10.filter,windowSize,...))%>%
dplyr::mutate(M.M21.T8.tricube=tricubeStat(POS=POS,Stat=M.M21.T8.filter,windowSize,...))%>%
dplyr::mutate(N.M1.tricube=tricubeStat(POS=POS,Stat=N.M1.filter,windowSize,...))%>%
dplyr::mutate(N.M1.T0.2.tricube=tricubeStat(POS=POS,Stat=N.M1.T0.2.filter,windowSize,...))%>%
dplyr::mutate(N.M11.T10.tricube=tricubeStat(POS=POS,Stat=N.M11.T10.filter,windowSize,...))%>%
dplyr::mutate(N.M11.T5.tricube=tricubeStat(POS=POS,Stat=N.M11.T5.filter,windowSize,...))%>%
dplyr::mutate(N.M11.T8.tricube=tricubeStat(POS=POS,Stat=N.M11.T8.filter,windowSize,...))%>%
dplyr::mutate(N.M2.tricube=tricubeStat(POS=POS,Stat=N.M2.filter,windowSize,...))%>%
dplyr::mutate(N.M2.T0.2.tricube=tricubeStat(POS=POS,Stat=N.M2.T0.2.filter,windowSize,...))%>%
dplyr::mutate(N.M21.T10.tricube=tricubeStat(POS=POS,Stat=N.M21.T10.filter,windowSize,...))%>%
dplyr::mutate(N.M21.T2.tricube=tricubeStat(POS=POS,Stat=N.M21.T2.filter,windowSize,...))%>%
dplyr::mutate(N.M21.T5.tricube=tricubeStat(POS=POS,Stat=N.M21.T5.filter,windowSize,...))%>%
dplyr::mutate(N.M21.T8.tricube=tricubeStat(POS=POS,Stat=N.M21.T8.filter,windowSize,...))%>%
dplyr::mutate(N.M3.tricube=tricubeStat(POS=POS,Stat=N.M3.filter,windowSize,...))%>%
dplyr::mutate(N.M3.T0.2.tricube=tricubeStat(POS=POS,Stat=N.M3.T0.2.filter,windowSize,...))%>%
dplyr::mutate(N.M4.tricube=tricubeStat(POS=POS,Stat=N.M4.filter,windowSize,...))%>%
dplyr::mutate(N.M4.T0.2.tricube=tricubeStat(POS=POS,Stat=N.M4.T0.2.filter,windowSize,...))
return(as.data.frame(SNPset))};refFre.AD.BSA6.LC<-SAnalysis.1(refFre.AD.BSA6.LC)


> dim(refFre.AD.BSA6.LC)
[1] 12803   221

write.csv(refFre.AD.BSA6.LC, file = "BSA6.ALL.AF.filter.tricube_coregenome.csv",row.names=FALSE)


################################################################
################################################################
################################################################

refFre.AD.BSA6.LC <- read.csv(file = "BSA6.ALL.AF.filter.tricube_coregenome.csv",header = TRUE)
refFre.AD.BSA6.LC$CHROM <- gsub("\\_v3|Pf3D7\\_0|Pf3D7\\_","",refFre.AD.BSA6.LC$CHROM)
refFre.AD.BSA6.LC$CHROM <- as.numeric(as.character(refFre.AD.BSA6.LC$CHROM))

################################################################
# plot allele frequency : life cycle 

p0 <- ggplot(data=refFre.AD.BSA6.LC)+scale_x_continuous(breaks=seq(from=0,to=max(refFre.AD.BSA6.LC$POS),by=10^(floor(log10(max(refFre.AD.BSA6.LC$POS))))),labels=format_genomic())+ylim(0,1)+facet_grid(~CHROM,scales="free_x",space="free_x") +theme_classic()

p0 + geom_point(aes_string(x = "POS", y = "C37.M11.D21.tricube"),color = "grey",size=0.5) + geom_line(aes_string(x = "POS", y = "C39.M11.D21.tricube"),color = "black",size=0.5)


# plot allele frequency : Temp BSA
library(ggpubr)

Temp.M11 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "C37.M11.D3.tricube"),color = "#74C476",size=0.5)+ ## 37C.day3
		 geom_line(aes_string(x = "POS", y = "C37.M11.D10.tricube"),color = "#11723B",size=0.5)+ ## 37C.day10
		 geom_line(aes_string(x = "POS", y = "C39.M11.D3.tricube"),color = "#FFEDA0",size=0.5)+ # 39C.day3
		 geom_line(aes_string(x = "POS", y = "C39.M11.D10.tricube"),color = "#FC6D33",size=0.5) # 39C.day10
		 #geom_line(aes_string(x = "POS", y = "C37.M11.D21.tricube"),color = "red",size=0.5)+ ## 37C.day21 too skew
		 #geom_line(aes_string(x = "POS", y = "C39.M11.D21.tricube"),color = "red",size=0.5,linetype = "dotted") # 39C.day21 too skew
Temp.M21 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "C37.M21.D3.tricube"),color = "#74C476",size=0.5)+ ## 37C.day3
		 geom_line(aes_string(x = "POS", y = "C37.M21.D10.tricube"),color = "#11723B",size=0.5)+ ## 37C.day10
		 geom_line(aes_string(x = "POS", y = "C39.M21.D3.tricube"),color = "#FFEDA0",size=0.5)+ # 39C.day3
		 geom_line(aes_string(x = "POS", y = "C39.M21.D10.tricube"),color = "#FC6D33",size=0.5) # 39C.day10
		 #geom_line(aes_string(x = "POS", y = "C37.M21.D21.tricube"),color = "red",size=0.5,linetype = "dotted")+ ## 37C.day21 too skew
		 #geom_line(aes_string(x = "POS", y = "C39.M21.D21.tricube"),color = "red",size=0.5,linetype = "dotted") # 39C.day21 too skew
		 		 
pdf('Temp BSA.pdf', width=10, height=6)
ggarrange(Temp.M11, Temp.M21,
          labels = c("Temp.M11", "Temp.M21"),
          ncol = 1, nrow = 2)
dev.off()

# Texas initial allele frequency

Nutri.T0.2 <- p0 +ylab("Mal31 allele frequency")+
          geom_line(aes_string(x = "POS", y = "N.M1.T0.2.tricube"),color = "#A1D99B",size=0.5)+ ## 
		  geom_line(aes_string(x = "POS", y = "N.M2.T0.2.tricube"),color = "orange",size=0.5)+ ## 
		  geom_line(aes_string(x = "POS", y = "N.M3.T0.2.tricube"),color = "red",size=0.5)+ ## 
		  geom_line(aes_string(x = "POS", y = "N.M4.T0.2.tricube"),color = "black",size=0.5) ## 
		  
Nutri.T0.1 <- p0 +ylab("Mal31 allele frequency")+
          geom_line(aes_string(x = "POS", y = "N.M1.tricube"),color = "#A1D99B",size=0.5)+ ## 
		  geom_line(aes_string(x = "POS", y = "N.M2.tricube"),color = "orange",size=0.5)+ ## 
		  geom_line(aes_string(x = "POS", y = "N.M3.tricube"),color = "red",size=0.5)+ ## 
		  geom_line(aes_string(x = "POS", y = "N.M4.tricube"),color = "black",size=0.5) ## 
		  		  
pdf('Texas BSA T0.pdf', width=10, height=6)
ggarrange(Nutri.T0.1, Nutri.T0.2,
          labels = c("Nutri.T0.1", "Nutri.T0.2"),
          ncol = 1, nrow = 2)
dev.off()
		  

# plot allele frequency : Nutri BSA
Nutri.M11 <- p0 + ylab("Mal31 allele frequency")+
         # geom_line(aes_string(x = "POS", y = "N.M11.T2.tricube"),color = "#D6EFD0",size=0.5)+ ## N.day2 too low coverage
		 geom_line(aes_string(x = "POS", y = "N.M11.T5.tricube"),color = "#A1D99B",size=0.5)+ ## N.day5
         geom_line(aes_string(x = "POS", y = "N.M11.T8.tricube"),color = "#74C476",size=0.5)+ ## N.day8
		 geom_line(aes_string(x = "POS", y = "N.M11.T10.tricube"),color = "#319B50",size=0.5)+ ## N.day10	
         geom_line(aes_string(x = "POS", y = "D.M11.T2.tricube"),color = "#FFEDA0",size=0.5)+ ## D.day2 
		 #geom_line(aes_string(x = "POS", y = "D.M11.T5.tricube"),color = "#FED976",size=0.5)+ ## D.day5 too low coverage
         geom_line(aes_string(x = "POS", y = "D.M11.T8.tricube"),color = "#FD9F44",size=0.5)+ ## D.day8
		 geom_line(aes_string(x = "POS", y = "D.M11.T10.tricube"),color = "#E31A1C",size=0.5)+ ## D.day10			 
         geom_line(aes_string(x = "POS", y = "M.M11.T8.tricube"),color = "grey",size=0.5)+ ## M.day8
		 geom_line(aes_string(x = "POS", y = "M.M11.T10.tricube"),color = "black",size=0.5) ## M.day10		 
		 
         #geom_line(aes_string(x = "POS", y = "M.M11.T2.tricube"),color = "#FFEDA0",size=0.5)+ ## M.day2 too low coverage
		 #geom_line(aes_string(x = "POS", y = "M.M11.T5.tricube"),color = "#FED976",size=0.5)+ ## M.day5 too low coverage
         geom_line(aes_string(x = "POS", y = "M.M11.T8.tricube"),color = "#FD9F44",size=0.5)+ ## M.day8
		 geom_line(aes_string(x = "POS", y = "M.M11.T10.tricube"),color = "#E31A1C",size=0.5) ## M.day10
		 

Nutri.M21 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "N.M21.T2.tricube"),color = "#D6EFD0",size=0.5)+ ## N.day2 too low coverage
		 geom_line(aes_string(x = "POS", y = "N.M21.T5.tricube"),color = "#A1D99B",size=0.5)+ ## N.day5
         geom_line(aes_string(x = "POS", y = "N.M21.T8.tricube"),color = "#74C476",size=0.5)+ ## N.day8
		 geom_line(aes_string(x = "POS", y = "N.M21.T10.tricube"),color = "#319B50",size=0.5)+ ## N.day10	
         geom_line(aes_string(x = "POS", y = "D.M21.T2.tricube"),color = "#FFEDA0",size=0.5)+ ## D.day2 
		 geom_line(aes_string(x = "POS", y = "D.M21.T5.tricube"),color = "#FED976",size=0.5)+ ## D.day5 too low coverage
         geom_line(aes_string(x = "POS", y = "D.M21.T8.tricube"),color = "#FD9F44",size=0.5)+ ## D.day8
		 geom_line(aes_string(x = "POS", y = "D.M21.T10.tricube"),color = "#E31A1C",size=0.5)+ ## D.day10			 
         geom_line(aes_string(x = "POS", y = "M.M21.T8.tricube"),color = "grey",size=0.5)+ ## M.day8
		 geom_line(aes_string(x = "POS", y = "M.M21.T10.tricube"),color = "black",size=0.5) ## M.day10		 
		 
         #geom_line(aes_string(x = "POS", y = "M.M21.T2.tricube"),color = "#FFEDA0",size=0.5)+ ## M.day2 too low coverage
		 #geom_line(aes_string(x = "POS", y = "M.M21.T5.tricube"),color = "#FED976",size=0.5)+ ## M.day5 too low coverage
         geom_line(aes_string(x = "POS", y = "M.M21.T8.tricube"),color = "#FD9F44",size=0.5)+ ## M.day8
		 geom_line(aes_string(x = "POS", y = "M.M21.T10.tricube"),color = "#E31A1C",size=0.5) ## M.day10

pdf('Nutri BSA.pdf', width=10, height=6)
ggarrange(Nutri.M11, Nutri.M21,
          labels = c("Nutri.M11", "Nutri.M21"),
          ncol = 1, nrow = 2)
dev.off()


# plot allele frequency : life cycle BSA
LC.1.albumax <- p0 + ylab("Mal31 allele frequency")+
         # geom_line(aes_string(x = "POS", y = "FG.BSA.0007.tricube"),color = "#FFEDA0",size=0.5)+ ## SGSpz mess
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0016.tricube"),color = "#FEE38B",size=0.5)+ ## Blood(in vivo)
         geom_line(aes_string(x = "POS", y = "FG.BSA.0017.tricube"),color = "#FEC561",size=0.5)+ ## Blood(in vivo)
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0043.tricube"),color = "#FEB24C",size=0.5)+ ## Blood day 7 Albumax-red
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0049.tricube"),color = "#FD9F44",size=0.5)+ ## Blood day 7 Albumax
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0067.tricube"),color = "#FD8D3C",size=0.5)+ ## Blood day 12 Albumax
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0115.tricube"),color = "#FC6D33",size=0.5)+ ## Blood day 18 Albumax
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0163.tricube"),color = "#FC4E2A",size=0.5)+ ## Blood day 24 Albumax
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0211.tricube"),color = "#EF3322",size=0.5) + ## Blood day 30 Albumax
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0259.tricube"),color = "#E31A1C",size=0.5) ## Blood day 36 Albumax

		 
LC.1.albumax_serum <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BSA.0067.tricube"),color = "#FD8D3C",size=0.5)+ ## Blood day 12 Albumax
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0115.tricube"),color = "#FC6D33",size=0.5)+ ## Blood day 18 Albumax
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0163.tricube"),color = "#FC4E2A",size=0.5)+ ## Blood day 24 Albumax
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0211.tricube"),color = "#EF3322",size=0.5) + ## Blood day 30 Albumax
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0259.tricube"),color = "#E31A1C",size=0.5) + ## Blood day 36 Albumax
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0075.tricube"),color = "#A1D99B",size=0.5)+ ## Blood day 12 Serum-green
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0123.tricube"),color = "#74C476",size=0.5)+ ## Blood day 18 Serum
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0171.tricube"),color = "#41AB5D",size=0.5)+ ## Blood day 24 Serum
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0219.tricube"),color = "#238B45",size=0.5)+ ## Blood day 30 Serum
		 geom_line(aes_string(x = "POS", y = "FG.BSA.0267.tricube"),color = "#005A32",size=0.5) ## Blood day 36 Serum
		 
pdf('lifecycle.serum.albumax BSA.pdf', width=10, height=6)
ggarrange(LC.1.albumax, LC.1.albumax_serum,
          labels = c("lifecycle-albuamx", "albumax vs serum"),
          ncol = 1, nrow = 2)
dev.off()


# plot allele frequency : drug BSA
##########################################
Drug.T0 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0196.tricube"),color = "#FFEDA0",size=0.5)+ ## 
		 geom_line(aes_string(x = "POS", y = "FG.BC.0197.tricube"),color = "#FEE38B",size=0.5)+ ## 
         geom_line(aes_string(x = "POS", y = "FG.BC.0198.tricube"),color = "#FEC561",size=0.5)+ ## 
		 geom_line(aes_string(x = "POS", y = "FG.BC.0199.tricube"),color = "#FEB24C",size=0.5)+ ## 
		 geom_line(aes_string(x = "POS", y = "FG.BC.0200.tricube"),color = "#FD9F44",size=0.5)+ ## 
		 geom_line(aes_string(x = "POS", y = "FG.BC.0201.tricube"),color = "#FD8D3C",size=0.5)+ ## 
		 geom_line(aes_string(x = "POS", y = "FG.BC.0202.tricube"),color = "#FC6D33",size=0.5)+ ## 
		 geom_line(aes_string(x = "POS", y = "FG.BC.0203.tricube"),color = "#FC4E2A",size=0.5)
		 
Drug.CQ.M3 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0196.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0258.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0261.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0264.tricube"),color = "red",size=0.5)+ ## 100nM day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0267.tricube"),color = "black",size=0.5) ## 250nM day5


Drug.CQ.M4 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0197.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0270.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0273.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0276.tricube"),color = "red",size=0.5)+ ## 100nM day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0279.tricube"),color = "black",size=0.5) ## 250nM day5


pdf('CQ BSA.pdf', width=10, height=6)
ggarrange(Drug.CQ.M3, Drug.CQ.M4,
          labels = c("CQ.M3", "CQ.M4"),
          ncol = 1, nrow = 2)
dev.off()

##########################################
DHA.M3 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0198.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0204.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0207.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0210.tricube"),color = "red",size=0.5) ## 100nM day5


DHA.M4 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0199.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0213.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0216.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0219.tricube"),color = "red",size=0.5)## 100nM day5



d1.2 <-importFromGATK(file = "BSA6.1.multi.MAL31xKH004.SNP.filter.table",highBulk = "FG.BC.0204",lowBulk = "FG.BC.0210")
d1.2.filt <-filterSNPs(SNPset = d1.2,minTotalDepth = 60,maxTotalDepth = 5000,minSampleDepth = 30,minGQ = 10)
d1.2.filt <- runQTLseqAnalysis(SNPset = d1.2.filt,windowSize = 1e5,popStruc = "F2",bulkSize = c(600, 600),replications = 10000,intervals = c(95, 99))

d1.2.filt <- runGprimeAnalysis(d1.2.filt,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)

plotQTLStats(SNPset = d1.2.filt, var = "Gprime", plotThreshold = FALSE, q = 0.05, line=TRUE,color = "orange",size= 1)+theme_classic()

write.csv(d1.2.filt, file = "DHA BSA_M3.csv",row.names=FALSE)




# for Manoj grant
setwd("C:/Users/xli.TXBIOMED/Dropbox (TX Biomed)/Emily/P01/6.BSA6.Mal31xKH004/3.analysis/BSA"): office desktop qtlseqr not working 
### DHA
QTL <- read.csv(file = "DHA BSA.csv",header = TRUE)
QTL$CHROM <- gsub("\\_v3|Pf3D7\\_0|Pf3D7\\_","",QTL$CHROM)
QTL$CHROM <- as.numeric(as.character(QTL$CHROM))
QTL.filt <- QTL[complete.cases(QTL), ]

p0 <- ggplot(data=QTL.filt)+scale_x_continuous(breaks=seq(from=0,to=max(QTL$POS),by=10^(floor(log10(max(QTL$POS))))),labels=format_genomic())+facet_grid(~CHROM,scales="free_x",space="free_x") +theme_classic()

p.DHA.Gprime.Manoj <- p0 + ylab(" G'(DHA) ")+ xlab("Genomic Position (Mb)") + 
     geom_hline(yintercept = 20,color = "red",size = 0.5, linetype = "dashed") +
     geom_line(aes_string(x = "POS", y = "M3.Gprime"),color = "black",size=1) +
     geom_line(aes_string(x = "POS", y = "M4.Gprime"),color = "orange",size=1)
	 

pdf('DHA Manoj grant.pdf', width=10, height=3)
ggarrange(p.DHA.Gprime.Manoj,
          labels = c(""),
          ncol = 1, nrow = 1)
dev.off()
	 
	 
### CQ

p.CQ.Gprime.Manoj <- p0 + ylab(" G'(CQ) ")+ xlab("Genomic Position (Mb)") + 
     geom_hline(yintercept = 20,color = "red",size = 0.5, linetype = "dashed") +
     geom_line(aes_string(x = "POS", y = "M3.Gprime"),color = "black",size=1) +
     geom_line(aes_string(x = "POS", y = "M4.Gprime"),color = "orange",size=1)
	 

pdf('CQ Manoj grant.pdf', width=10, height=3)
ggarrange(p.CQ.Gprime.Manoj,
          labels = c(""),
          ncol = 1, nrow = 1)
dev.off()
	 
	 




pdf('DHA high dose BSA.pdf', width=10, height=6)
ggarrange(DHA.M3, DHA.M4,
          labels = c("DHA.M3", "DHA.M4"),
          ncol = 1, nrow = 2)
dev.off()

DHA.M3.day8 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0198.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0318.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0321.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0324.tricube"),color = "red",size=0.5) ## 100nM day5


DHA.M4.day8 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0199.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0327.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 # geom_line(aes_string(x = "POS", y = "FG.BC.0330.tricube"),color = "orange",size=0.5)+ ## 50nM day5 # too low coverage
         geom_line(aes_string(x = "POS", y = "FG.BC.0333.tricube"),color = "red",size=0.5)## 100nM day5

pdf('DHA high dose BSA - day 8.pdf', width=10, height=6)
ggarrange(DHA.M3.day8, DHA.M4.day8 ,
          labels = c("DHA.M3.d8", "DHA.M4.d8"),
          ncol = 1, nrow = 2)
dev.off()

##########################################
DHA.P18.M3 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0200.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0282.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0285.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0288.tricube"),color = "red",size=0.5)## 100nM day5
		 
DHA.P36.M3 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0202.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0300.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0303.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0306.tricube"),color = "red",size=0.5)## 100nM day5

pdf('DHA high dose BSA-synched.pdf', width=10, height=6)
ggarrange(DHA.P18.M3, DHA.P36.M3,
          labels = c("DHA.P18.M3", "DHA.P36.M3"),
          ncol = 1, nrow = 2)
dev.off()


##########################################
PPQ.M3 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0196.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0222.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0225.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0228.tricube"),color = "red",size=0.5) ## 100nM day5

PPQ.M3.Day8 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0196.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0336.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0339.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0342.tricube"),color = "red",size=0.5) ## 100nM day5


PPQ.M4 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0197.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0231.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0234.tricube"),color = "orange",size=0.5) ## 50nM day5
         #geom_line(aes_string(x = "POS", y = "FG.BC.0237.tricube"),color = "red",size=0.5) ## 100nM day5  # too low coverage


pdf('PPQ BSA.pdf', width=10, height=6)
ggarrange(PPQ.M3, PPQ.M4,
          labels = c("PPQ.M3", "PPQ.M4"),
          ncol = 1, nrow = 2)
dev.off()



##########################################
MB.M3 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0196.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0240.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0243.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0246.tricube"),color = "red",size=0.5) ## 100nM day5

MB.M4 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0197.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0249.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0252.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0255.tricube"),color = "red",size=0.5) ## 100nM day5


MB.M3.DAY8 <- p0 + ylab("Mal31 allele frequency")+
         geom_line(aes_string(x = "POS", y = "FG.BC.0198.tricube"),color = "green",size=1)+ ## before 
         geom_line(aes_string(x = "POS", y = "FG.BC.0354.tricube"),color = "#FFEDA0",size=0.5)+ ## control day5
		 geom_line(aes_string(x = "POS", y = "FG.BC.0357.tricube"),color = "orange",size=0.5)+ ## 50nM day5
         geom_line(aes_string(x = "POS", y = "FG.BC.0360.tricube"),color = "red",size=0.5) ## 100nM day5  # too low coverage


pdf('MB BSA.pdf', width=10, height=6)
ggarrange(MB.M3, MB.M4,
          labels = c("MB.M3", "MB.M4"),
          ncol = 1, nrow = 2)
dev.off()


