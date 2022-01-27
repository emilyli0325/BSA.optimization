
#### deltaSNP and Gprime
dfOo <-importFromGATK(file = "BSA.7.2.NF54xNHP1337.Drug.SNP.filter.table",highBulk = "BE.BC.0027",lowBulk = "BE.BC.0033")
df_filtOo <-filterSNPs(SNPset = dfOo,minTotalDepth = 30,maxTotalDepth = 10000,minSampleDepth = 10)
df_filtOo <- runQTLseqAnalysis(SNPset = df_filtOo,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filtOo <- runGprimeAnalysis(df_filtOo,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
plotQTLStats(SNPset = df_filtOo, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
plotQTLStats(SNPset = df_filtOo, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()+theme(strip.background = element_blank(),strip.text.x = element_blank())



BSA <- read.table("BSA.7.2.NF54xNHP1337.Drug.SNP.filter.table", sep = '\t',header = TRUE)


> dim(BSA)
[1] 16819    84


AD <- BSA[,seq(5,84,4)]
DP <- BSA[,seq(6,84,4)]
ref.DP <- function(X){as.numeric(strsplit(as.character(X),",")[[1]])[1]}
dim(AD)
[1] 16819    20
refFre.AD <- matrix(ncol=20,nrow=16819)
for(j in 1:20){
	refFre.AD[,j]<-sapply(AD[,j],ref.DP,simplify="array")/as.numeric(DP[,j])
	}
refFre.AD[DP<20]<-NA
colnames(refFre.AD)<-colnames(AD)	
refFre.AD <- cbind(BSA[,1:4],refFre.AD)
colnames(refFre.AD) <- gsub(".AD", "", colnames(refFre.AD))

### remove outliers ###

outliersMAD <- function(data, MADCutOff = 2.5, replace = NA, values = FALSE, bConstant = 1.4826, digits = 2) {
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


							   
## filter
AB.BSA.220.NF54.filter <- outlierByMAD(refFre.AD$AB.BSA.220.NF54,50)
BE.BC.0001.filter <- outlierByMAD(refFre.AD$BE.BC.0001,50)
BE.BC.0002.filter <- outlierByMAD(refFre.AD$BE.BC.0002,50)
BE.BC.0003.filter <- outlierByMAD(refFre.AD$BE.BC.0003,50)
BE.BC.0004.filter <- outlierByMAD(refFre.AD$BE.BC.0004,50)
BE.BC.0005.filter <- outlierByMAD(refFre.AD$BE.BC.0005,50)
BE.BC.0006.filter <- outlierByMAD(refFre.AD$BE.BC.0006,50)
BE.BC.0007.filter <- outlierByMAD(refFre.AD$BE.BC.0007,50)
BE.BC.0008.filter <- outlierByMAD(refFre.AD$BE.BC.0008,50)
BE.BC.0009.filter <- outlierByMAD(refFre.AD$BE.BC.0009,50)
BE.BC.0012.filter <- outlierByMAD(refFre.AD$BE.BC.0012,50)
BE.BC.0015.filter <- outlierByMAD(refFre.AD$BE.BC.0015,50)
BE.BC.0018.filter <- outlierByMAD(refFre.AD$BE.BC.0018,50)
BE.BC.0021.filter <- outlierByMAD(refFre.AD$BE.BC.0021,50)
BE.BC.0024.filter <- outlierByMAD(refFre.AD$BE.BC.0024,50)
BE.BC.0027.filter <- outlierByMAD(refFre.AD$BE.BC.0027,50)
BE.BC.0030.filter <- outlierByMAD(refFre.AD$BE.BC.0030,50)
BE.BC.0033.filter <- outlierByMAD(refFre.AD$BE.BC.0033,50)
BE.BC.0036.filter <- outlierByMAD(refFre.AD$BE.BC.0036,50)
NHP1337.12C.filter <- outlierByMAD(refFre.AD$NHP1337.12C,50)


refFre.AD.filter <- cbind(refFre.AD[,1:4], AB.BSA.220.NF54.filter,BE.BC.0001.filter,BE.BC.0002.filter,BE.BC.0003.filter,BE.BC.0004.filter,BE.BC.0005.filter,BE.BC.0006.filter,BE.BC.0007.filter,BE.BC.0008.filter,BE.BC.0009.filter,BE.BC.0012.filter,BE.BC.0015.filter,BE.BC.0018.filter,BE.BC.0021.filter,BE.BC.0024.filter,BE.BC.0027.filter,BE.BC.0030.filter,BE.BC.0033.filter,BE.BC.0036.filter,NHP1337.12C.filter)



tricubeStat <- function(POS, Stat, windowSize = 1e5, ...)
{
    if (windowSize <= 0)
        stop("A positive smoothing window is required")
    stats::predict(locfit::locfit(Stat ~ locfit::lp(POS, h = windowSize, deg = 0), ...), POS)
}

SAnalysis.1 <-
    function(SNPset,
        windowSize = 1e5,
        ...)
    {
 
        SNPset <- SNPset %>%
            dplyr::group_by(CHROM) %>%
            dplyr::mutate(AB.BSA.220.NF54.tricube = tricubeStat(POS = POS, Stat = AB.BSA.220.NF54.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0001.tricube = tricubeStat(POS = POS, Stat = BE.BC.0001.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0002.tricube = tricubeStat(POS = POS, Stat = BE.BC.0002.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0003.tricube = tricubeStat(POS = POS, Stat = BE.BC.0003.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0004.tricube = tricubeStat(POS = POS, Stat = BE.BC.0004.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0005.tricube = tricubeStat(POS = POS, Stat = BE.BC.0005.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0006.tricube = tricubeStat(POS = POS, Stat = BE.BC.0006.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0007.tricube = tricubeStat(POS = POS, Stat = BE.BC.0007.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0008.tricube = tricubeStat(POS = POS, Stat = BE.BC.0008.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0009.tricube = tricubeStat(POS = POS, Stat = BE.BC.0009.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0012.tricube = tricubeStat(POS = POS, Stat = BE.BC.0012.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0015.tricube = tricubeStat(POS = POS, Stat = BE.BC.0015.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0018.tricube = tricubeStat(POS = POS, Stat = BE.BC.0018.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0021.tricube = tricubeStat(POS = POS, Stat = BE.BC.0021.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0024.tricube = tricubeStat(POS = POS, Stat = BE.BC.0024.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0027.tricube = tricubeStat(POS = POS, Stat = BE.BC.0027.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0030.tricube = tricubeStat(POS = POS, Stat = BE.BC.0030.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0033.tricube = tricubeStat(POS = POS, Stat = BE.BC.0033.filter, windowSize, ...))  %>%
            dplyr::mutate(BE.BC.0036.tricube = tricubeStat(POS = POS, Stat = BE.BC.0036.filter, windowSize, ...))  %>%
            dplyr::mutate(NHP1337.12C.tricube = tricubeStat(POS = POS, Stat = NHP1337.12C.filter, windowSize, ...))
		return(as.data.frame(SNPset))	
    }

refFre.AD.tri <- SAnalysis.1(refFre.AD.filter)		


	
p0 <- ggplot(data = refFre.AD.tri) +scale_x_continuous(breaks = seq(from = 0,to = max(refFre.AD.tri$POS), by = 10^(floor(log10(max(refFre.AD.tri$POS))))), labels=format_genomic())+ facet_grid(~ CHROM, scales = "free_x", space = "free_x")+theme_classic()


p01 <- p0 + ylab("NF54 allele frequency") +
            geom_point(aes_string(x = "POS", y = "BE.BC.0001.filter"),color = "grey",size=0.5)+ 
            geom_point(aes_string(x = "POS", y = "BE.BC.0005.filter"),color = "orange",size=0.5)+
            geom_line(aes_string(x = "POS", y = "BE.BC.0001.tricube"),color = "black",size=1)+
		    geom_line(aes_string(x = "POS", y = "BE.BC.0005.tricube"),color = "red",size=1)+
		    geom_hline(yintercept = 0.5,color = "black",linetype="dashed",alpha = 1,size=1)
p02 <- p0 + ylab("NF54 allele frequency") +
            geom_point(aes_string(x = "POS", y = "BE.BC.0002.filter"),color = "grey",size=0.5)+ 
            geom_point(aes_string(x = "POS", y = "BE.BC.0006.filter"),color = "orange",size=0.5)+
            geom_line(aes_string(x = "POS", y = "BE.BC.0002.tricube"),color = "black",size=1)+
		    geom_line(aes_string(x = "POS", y = "BE.BC.0006.tricube"),color = "red",size=1)+
		    geom_hline(yintercept = 0.5,color = "black",linetype="dashed",alpha = 1,size=1)
p03 <- p0 + ylab("NF54 allele frequency") +
            geom_point(aes_string(x = "POS", y = "BE.BC.0003.filter"),color = "grey",size=0.5)+ 
            geom_point(aes_string(x = "POS", y = "BE.BC.0007.filter"),color = "orange",size=0.5)+
            geom_line(aes_string(x = "POS", y = "BE.BC.0003.tricube"),color = "black",size=1)+
		    geom_line(aes_string(x = "POS", y = "BE.BC.0007.tricube"),color = "red",size=1)+
		    geom_hline(yintercept = 0.5,color = "black",linetype="dashed",alpha = 1,size=1)
p04 <- p0 + ylab("NF54 allele frequency") +
            geom_point(aes_string(x = "POS", y = "BE.BC.0004.filter"),color = "grey",size=0.5)+ 
            geom_point(aes_string(x = "POS", y = "BE.BC.0008.filter"),color = "orange",size=0.5)+
            geom_line(aes_string(x = "POS", y = "BE.BC.0004.tricube"),color = "black",size=1)+
		    geom_line(aes_string(x = "POS", y = "BE.BC.0008.tricube"),color = "red",size=1)+
		    geom_hline(yintercept = 0.5,color = "black",linetype="dashed",alpha = 1,size=1)
		   

pdf('NF54xNHP1337_T0.pdf', width=10, height=10)
ggarrange(p01, p02, p03, p04,labels = c("A","B","C","D"), ncol = 1, nrow = 4)
dev.off()

## DHA BSA

p0 <- ggplot(data = refFre.AD.tri) +scale_x_continuous(breaks = seq(from = 0,to = max(refFre.AD.tri$POS), by = 10^(floor(log10(max(refFre.AD.tri$POS))))), labels=format_genomic())+ facet_grid(~ CHROM, scales = "free_x", space = "free_x")+theme_classic()


DHA.P0 <- p0 + ylab("NF54 allele frequency") +
            geom_line(aes_string(x = "POS", y = "BE.BC.0001.tricube"),color = "black",size=1)+
            geom_line(aes_string(x = "POS", y = "BE.BC.0009.tricube"),color = "green",size=1)+
		    geom_line(aes_string(x = "POS", y = "BE.BC.0012.tricube"),color = "orange",size=1)+
		    geom_line(aes_string(x = "POS", y = "BE.BC.0015.tricube"),color = "red",size=1)+
		    geom_hline(yintercept = 0.5,color = "black",linetype="dashed",alpha = 1,size=1)
			
DHA.unsync <- p0 + ylab("NF54 allele frequency") +
            geom_line(aes_string(x = "POS", y = "BE.BC.0005.tricube"),color = "black",size=1)+
            geom_line(aes_string(x = "POS", y = "BE.BC.0018.tricube"),color = "green",size=1)+
		    geom_line(aes_string(x = "POS", y = "BE.BC.0021.tricube"),color = "orange",size=1)+
			geom_line(aes_string(x = "POS", y = "BE.BC.0024.tricube"),color = "red",size=1)+
		    geom_hline(yintercept = 0.5,color = "black",linetype="dashed",alpha = 1,size=1)
			
			
pdf('NF54xNHP1337_DHA.BSA.pdf', width=10, height=5)
ggarrange(DHA.P0, DHA.unsync,labels = c("DHA.P0","DHA.unsync"), ncol = 1, nrow = 2)
dev.off()


## CQ BSA

p0 <- ggplot(data = refFre.AD.tri) +scale_x_continuous(breaks = seq(from = 0,to = max(refFre.AD.tri$POS), by = 10^(floor(log10(max(refFre.AD.tri$POS))))), labels=format_genomic())+ facet_grid(~ CHROM, scales = "free_x", space = "free_x")+theme_classic()


CQ.unsync <- p0 + ylab("NF54 allele frequency") +
            geom_line(aes_string(x = "POS", y = "BE.BC.0001.tricube"),color = "black",size=1)+
            geom_line(aes_string(x = "POS", y = "BE.BC.0027.tricube"),color = "green",size=1)+
		    geom_line(aes_string(x = "POS", y = "BE.BC.0030.tricube"),color = "yellow",size=1)+
		    geom_line(aes_string(x = "POS", y = "BE.BC.0033.tricube"),color = "orange",size=1)+
		    geom_line(aes_string(x = "POS", y = "BE.BC.0036.tricube"),color = "red",size=1)+			
		    geom_hline(yintercept = 0.5,color = "black",linetype="dashed",alpha = 1,size=1)

			
pdf('NF54xNHP1337_cq.BSA.pdf', width=10, height=2.5)
ggarrange(CQ.unsync,labels = c("CQ.unsync"), ncol = 1, nrow = 1)
dev.off()


###### Gprime

dfOo <-importFromGATK(file = "BSA.7.2.NF54xNHP1337.Drug.SNP.filter.table",highBulk = "BE.BC.0030",lowBulk = "BE.BC.0036")
df_filtOo <-filterSNPs(SNPset = dfOo,minTotalDepth = 60,maxTotalDepth = 1000,minSampleDepth = 30)
df_filtOo <- runQTLseqAnalysis(SNPset = df_filtOo,windowSize = 1e5,popStruc = "F2",bulkSize = c(50, 50),replications = 10000,intervals = c(95, 99))
df_filtOo <- runGprimeAnalysis(df_filtOo,windowSize = 1e5,outlierFilter = "deltaSNP",filterThreshold = 0.05)
p1 <- plotQTLStats(SNPset = df_filtOo, var = "Gprime", plotThreshold = TRUE, q = 0.01, line=TRUE,color = "orange",size= 1)+theme_classic()
p2 <- plotQTLStats(SNPset = df_filtOo, var = "deltaSNP", line=TRUE,color = "orange",size=1, plotIntervals = TRUE)+theme_classic()+theme(strip.background = element_blank(),strip.text.x = element_blank())
p1
p2


write.csv(df_filtOo, "CQ.250nM.csv")






