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
outputSuffix <- c("_CumuInc.pdf")

geneStrings <- c("HLX", "KALRN", "CLSTN2_2", "RAB2B1M", "RAB2B2", "LOC107986324", "SEMA6A", "TNKS_2", "GRIK4", "CORO1C", "CALM1", "TENM1")
labelStrings <- c("HLX", "KALRN", "CLSTN2", "RAP2B (SNP01)", "RAP2B (SNP02)", "LOC107986324", "SEMA6A", "TNKS", "GRIK4", "CORO1C", "CALM1", "TENM1")
wtAllele <- c("GG", "GG", "TT", "GG", "CC", "CC", "TT", "AA", "AA", "TT", "GG", "GG")
mtAllele <- c("GA", "GC", "TC", "GA", "CG", "CT", "TC", "AG", "AG", "TC", "GA", "GA")
hmAllele <- c("AA", "CC", "CC", "AA", "GG", "TT", "CC", "GG", "GG", "CC", "AA", "AA")
for (i in 1:length(geneStrings)) {
  goi <- geneStrings[i]
  loi <- labelStrings[i]
  # SUBSET ONLY AFFECTED INDIVIDUALS #
  inputData_temp <- subset(inputData, !is.na(inputData[[goi]]) & inputData$LRRK2 > 0)
  inputData_temp[[goi]][inputData_temp[[goi]] == 0] = "Non-Carriers" #wtAllele[i]
  inputData_temp[[goi]][inputData_temp[[goi]] == 1] = "Carriers" #paste(mtAllele[i], hmAllele[i], sep=" + ")
  inputData_temp[[goi]][inputData_temp[[goi]] == 2] = "Carriers" #paste(mtAllele[i], hmAllele[i], sep=" + ")
  
  #compriskEst <- tidycmprsk::crr(Surv(Time, as.factor(Affected)) ~ CORO1C, data=inputData_CORO1C)
  callStr <- paste("Surv(Time, as.factor(Affected)) ~ ", goi, sep="")
  compriskFit <- tidycmprsk::cuminc(as.formula(callStr), data=inputData_temp)
  print(compriskFit)
  
  outputName <- paste(outputDirectory, outputPrefix, loi, outputSuffix, sep="")
  
  compriskFit %>%
    ggcuminc() + 
    add_confidence_interval() +
    add_risktable(risktable_stats = "{n.risk} ({cum.event})",
                  stats_label = list(n.risk = "Number at Risk", cum.event = "Number Affected")) + 
#                  ylab = c("1", "2")) + 
    scale_color_manual(name="Genotype", values=c("#696969", "black"), guide="legend", labels=c(paste(mtAllele[i], hmAllele[i], sep=" + "), wtAllele[i])) +
    scale_fill_manual(name="Genotype", values=c("#696969", "black"), guide="legend", labels=c(paste(mtAllele[i], hmAllele[i], sep=" + "), wtAllele[i])) +
    scale_x_continuous(limits=c(0, 100)) + 
    labs(title=paste(loi, " on Risk", sep=""), x="Time (Years)", y="Cumulative Incidence") + 
    theme(legend.position="bottom") +
    theme_minimal()
  #dev.off()
  ggsave(outputName, width=14, height=8.5)
}


