#
#
#

library(cmprsk)
library(readxl)
library(ggsurvfit)
library(stringr)

strictNoMutCohortPath <- "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\inputData_NoMutationStrictCohort_Updated-08-25-23.xlsx"
liberalNoMutCohortPath <- "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\inputData_NoMutationLiberalCohort_Updated-08-25-23.xlsx"

inputData <- read_excel(liberalNoMutCohortPath)

outputDirectory <- c("C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\Figures\\")
outputPrefix <- c("TunisiaLRRK2Modifiers-")
outputSuffix <- c("_AAO-Boxplots.pdf")

geneStrings <- c("HLX", "KALRN", "CLSTN2_2", "RAB2B1M", "RAB2B2", "LOC107986324", "SEMA6A", "TNKS_2", "GRIK4", "CORO1C", "CALM1", "TENM1")
labelStrings <- c("HLX", "KALRN", "CLSTN2", "RAP2B (SNP01)", "RAP2B (SNP02)", "LOC107986324", "SEMA6A", "TNKS", "GRIK4", "CORO1C", "CALM1", "TENM1")
wtAllele <- c("GG", "GG", "TT", "GG", "CC", "CC", "TT", "AA", "AA", "TT", "GG", "GG")
mtAllele <- c("GA", "GC", "TC", "GA", "CG", "CT", "TC", "AG", "AG", "TC", "GA", "GA")
hmAllele <- c("AA", "CC", "CC", "AA", "GG", "TT", "CC", "GG", "GG", "CC", "AA", "AA")

fetchSEM <- function(x){
  v <- c(min(x), mean(x) - (sd(x)/sqrt(length(x))), mean(x), mean(x) + (sd(x)/sqrt(length(x))), max(x))
  names(v) <- c("ymin", "lower", "middle", "upper", "ymax")
  return(v)
}

for (i in 1:length(geneStrings)) {
  goi <- geneStrings[i]
  loi <- labelStrings[i]
  mutIDX <- which(!is.na(inputData[[goi]]) & inputData[[goi]] > 0 & inputData$LRRK2 > 0)
  norIDX <- which(!is.na(inputData[[goi]]) & inputData[[goi]] == 0 & inputData$LRRK2 > 0)
  
  mutGSAffected <- inputData$Time[mutIDX]
  norGSAffected <- inputData$Time[norIDX]
  dot <- data.frame(
    data = c(mutGSAffected, norGSAffected), 
    group = c(rep("MT", length(mutGSAffected)), rep("WT", length(norGSAffected)))
  )

  outputName <- paste(outputDirectory, outputPrefix, loi, outputSuffix, sep="")
  
  dot %>%
    ggplot( aes(x=group, y=data, fill=group)) +
      stat_summary(fun.data=fetchSEM, geom="boxplot", color="black", fill=c("#C0C0C0", "#696969")) +
      stat_boxplot(geom="errorbar", width=0.6, coef=NULL) + 
      scale_x_discrete(labels = c(paste(mtAllele[i], hmAllele[i], sep=" + "), wtAllele[i])) + 
      scale_fill_discrete(labels = c(paste(mtAllele[i], hmAllele[i], sep=" + "), wtAllele[i])) +
      #geom_boxplot(color="black", fill=c("#C0C0C0", "#696969"), alpha=0.6, outlier.shape=NA) + 
      geom_jitter(color="black", size=1.6, alpha=1, width=0.3) + 
      ggtitle(paste(loi, " Genotype on AAO in Affected p.G2019S Carriers", sep="")) +
      xlab("Genotype") +
      ylab("Age (years)") +
      coord_cartesian(ylim = c(0, 100)) + 
      theme_bw()
    
  ggsave(outputName, width=14, height=8.5)
  
  resu <- t.test(data ~ group, dot)
  print(resu)
  mn1 <- mean(dot$data[dot$group == paste(mtAllele[i], hmAllele[i], sep=" + ") & !is.na(dot$data)])
  sd1 <- sd(dot$data[dot$group == paste(mtAllele[i], hmAllele[i], sep=" + ") & !is.na(dot$data)])
  ld1 <- length((dot$data[dot$group == paste(mtAllele[i], hmAllele[i], sep=" + ") & !is.na(dot$data)]))
  sem1 <- sd1/sqrt(ld1)
  mn2 <- mean(dot$data[dot$group == wtAllele[i] & !is.na(dot$data)])
  sd2 <- sd(dot$data[dot$group == wtAllele[i] & !is.na(dot$data)])
  ld2 <- length(dot$data[dot$group == wtAllele[i] & !is.na(dot$data)])
  sem2 <- sd2/sqrt(ld2)
  print(paste(loi, ": MT-", sem1, "(", ld1, "), WT-", sem2, "(", ld2, ")", sep=""))
}
