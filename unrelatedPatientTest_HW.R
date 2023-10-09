if (!"genetics" %in% installed.packages()) {
  install.packages("genetics")
}
library(genetics)
if (!"HardyWeinberg" %in% installed.packages()) {
  install.packages("HardyWeinberg")
}
library(HardyWeinberg)

# ---------------------------------------- #
# USE THIS SECTION TO ASSIGN ALL CONSTANTS #
# ---------------------------------------- #
mainLabels <- c("HLX,DUSP10", "KALRN", "CLSTN2", "RAP2B,ARHGEF26", "RAP2B,ARHGEF26", "LOC107986324", "SEMA6A", "TNKS", "NR6A1", "GRIK4", "LRRK2", "CORO1C", "CALM1,TTC7B", "TENM1")
mainRSIDs <- c("rs141686162", "rs145611031", "rs150382576", "rs59679443", "rs16846845", "rs12272007", "rs73781088", "rs28398294", "rs148922482", "rs28470321", "rs34637584", "rs77395454", "rs76788674", "rs185981774")

reduLabels <- c("", "KALRN", "CLSTN2", "RAP2B_A", "RAP2B_B", "LOC107986324", "SEMA6A", "TNKS", "NR6A1", "GRIK4", "LRRK2", "CORO1C", "CALM1", "TENM1")
reduRSIDs <- c("rs115340641", "rs114392312", "rs140825591", "rs73007728", "rs58923546", "rs6824505", "rs17509902", "rs28511573", "rs150792641", "rs10502241", "rs190882284", "rs3825252", "rs76039170", "rs185078407")

workspaceDir <- "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\"
#workspaceDir <- "C:\\Users\\HAL2\\Downloads\\"
agenaDataName <- "08.08.23_Tunisia_All_Info_SNP_Calls_DylanEdit.csv"
agenaDataDir <- paste(workspaceDir, agenaDataName, sep = "")

pedFileName <- "Tunisians linkage (format (UPN).PED"
pedFileDir <- "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\"
pedFileFull <- paste(pedFileDir, pedFileName, sep="")


# --------------------------------------------------------------- #
# THIS SECTION IS WHERE THE RAW DATA IS LOADED INTO OUR WORKSPACE #
# --------------------------------------------------------------- #
agenaData <- read.csv(agenaDataDir)
pedData <- read.table(pedFileFull, sep=" ", header=FALSE)


# -------------------------------------------------------------------------- #
# SUBSET SAMPLES FROM PROGENY LIST THAT HAVE A ZERO FOR BOTH MOM AND DAD IDS #
# -------------------------------------------------------------------------- #
potentialIDX <- which(pedData$V3 == 0 & pedData$V4 == 0 & pedData$V12 == '', arr.ind = TRUE)
potentialPats <- as.character(pedData$V9[potentialIDX])
agenaPats <- as.character(unique(agenaData$Individual.name))
safePats <- c()

for (i in potentialPats) {
  if (any(agenaPats == i)) {
    safePats <- c(safePats, i)
  }
}

mCount <- 0
fCount <- 0

for (i in 1:length(safePats)) {
  idxMoF <- which(pedData$V9 == safePats[i])
  if (pedData$V10[idxMoF] == 'F') {
    fCount <- fCount + 1
  } else if (pedData$V10[idxMoF] == 'M') {
    mCount <- mCount + 1
  }
}

numFemales <- length(which(agenaData$Gender == "F"))

for (i in 1:length(mainRSIDs)) {
  genoDataM <- c()
  genoDataR <- c()
  LDDataM <- c()
  LDDataR <- c()
  genoOIM <- c()
  genoOIR <- c()
  dataOIM <- subset(agenaData, agenaData$Assay == mainRSIDs[i])
  dataOIR <- subset(agenaData, agenaData$Assay == reduRSIDs[i])
  if (i == 14) {
    dataOIM <- subset(dataOIM, dataOIM$Gender == "F")
    dataOIR <- subset(dataOIR, dataOIR$Gender == "F")
  }
  for (j in safePats) {
    if (any(dataOIM$Individual.name == j)) {
      genoOIM <- as.character(subset(dataOIM$Replicates.Result, dataOIM$Individual.name == j))
      if (genoOIM != "mismatched") {
        if (nchar(genoOIM) == 1) {
          genoDataM <- append(genoDataM, paste(genoOIM, genoOIM, sep="/"))
        } else if (nchar(genoOIM) == 2) {
          genoDataM <- append(genoDataM, paste(substr(genoOIM, 1, 1), substr(genoOIM, 2, 2), sep="/"))
        } else {
          print(paste("POTENTIAL ERROR: geno longer than 2 in ", mainRSIDs[i], sep=""))
        }
      }
    }
    
    if (any(dataOIR$Individual.name == j)) {
      genoOIR <- as.character(subset(dataOIR$Replicates.Result, dataOIR$Individual.name == j))
      if (reduRSIDs[i] != "rs150792641") {
        if (genoOIR != "mismatched") {
          if (nchar(genoOIR) == 1) {
            genoDataR <- append(genoDataR, paste(genoOIR, genoOIR, sep="/"))
          } else if (nchar(genoOIR) == 2) {
            genoDataR <- append(genoDataR, paste(substr(genoOIR, 1, 1), substr(genoOIR, 2, 2), sep="/"))
          } else {
            print(paste("POTENTIAL ERROR: geno longer than 2 in ", reduRSIDs[i], sep=""))
          }
        }
      } else {
      
      }
    }
    
    if (length(genoOIM) > 0 & length(genoOIR) > 0) {
      if ((genoOIM != "mismatched" & genoOIR != "mismatched") & (any(dataOIM$Individual.name == j) & any(dataOIR$Individual.name == j))) {
        LDDataM <- append(LDDataM, genoDataM[length(genoDataM)])
        LDDataR <- append(LDDataR, genoDataR[length(genoDataR)])
      }
    }
  }
  genoOutM <- genotype(genoDataM)
  summary(genoOutM)
  HWE.chisq(genoOutM)
  genoOutR <- genotype(genoDataR)
  summary(genoOutR)
  HWE.chisq(genoOutR)
  
  LDDataM <- genotype(LDDataM)
  LDDataR <- genotype(LDDataR)
  LD(LDDataM, LDDataR)
  length(LDDataM)
}
