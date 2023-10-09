# CoxME test with kinship correction script
# For now: assuming Coro1c only modifier being run
# Model currently is: AOO = gender + G2019S_Geno + CORO1C_Geno + Affected_Status

if (!"coxme" %in% installed.packages()) {
  install.packages("coxme")
}
if (!"kinship2" %in% installed.packages()) {
  install.packages("kinship2")
}
if (!"survminer" %in% installed.packages()) {
  install.packages("survminer")
}
if (!"survival" %in% installed.packages()) {
  install.packages("survival")
}
if (!"ggplot2" %in% installed.packages()) {
  install.packages("ggplot2")
}
if (!"mstate" %in% installed.packages()) {
  install.packages("mstate")
}
if (!"cmprsk" %in% installed.packages()) {
  install.packages("cmprsk")
}
if (!"writexl" %in% installed.packages()) {
  install.packages("writexl")
}
if (!"readxl" %in% installed.packages()) {
  install.packages("readxl")
}

library(coxme)
library(kinship2)
library(survminer)
library(mstate)
library(dplyr)
library(tibble)
library(cmprsk)
library(writexl)
library(readxl)

assayGenes <- c("HLX", "KALRN", "CLSTN2", "CLSTN2_2", "LOC107986324", "SEMA6A", "TNKS", "TNKS_2", "NR6A1", "CALM1", "TENM1", "LRRK2", "CORO1C", "GRIK4", "RAB2B(1M)", "RAB2B(1R)", "RAB2B(2)")
assayLabels <- c("rs141686162", "rs145611031", "rs150382576", "rs140825591", "rs6824505", "rs73781088", "rs28398294", "rs28511573", "rs150792641", "rs76788674", "rs185981774", "rs34637584", "rs77395454", "rs28470321", "rs59679443", "rs73007728", "rs16846845")
assayWT <- c("G", "G", "G", "T", "C", "T", "A", "A", "AG", "G", "G", "G", "T", "A", "G", "T", "C")
assayHet <- c("AG", "GC", "AG", "CT", "CT", "CT", "GA", "AG", "AGAG.AG", "GA", "AG", "AG", "CT", "GA", "GA", "TC", "GC")
assayHom <- c("A", "C", "A", "C", "T", "C", "G", "G", "AGAG", "A", "A", "A", "C", "G", "A", "C", "G")
assayCounts <- data.frame(matrix(0, nrow=length(assayGenes), ncol=6))

workspaceDir <- "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\"
agenaDataName <- "08.09.23_Tunisia_All_Info_SNP_Calls_DylanEdit.csv"
agenaDataDir <- paste(workspaceDir, agenaDataName, sep = "")

pedFileName <- "Tunisians linkage (format (UPN).PED"
pedFileDir <- "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\"
pedFileFull <- paste(pedFileDir, pedFileName, sep="")

lrrk2HetTag <- "LRRK2 p.G2019S+/-"
lrrk2HomTag <- "LRRK2 p.G2019S+/+"

# --------------------------------------------------------------- #
# THIS SECTION IS WHERE THE RAW DATA IS LOADED INTO OUR WORKSPACE #
# --------------------------------------------------------------- #
agenaData <- read.csv(agenaDataDir)
pedData <- read.table(pedFileFull, sep=" ", header=FALSE)


# ---------------------------------------------------------------------------- #
# SUBSET SAMPLES FROM PROGENY LIST THAT HAVE AGENA RESULTS FOR SNP OF INTEREST #
# ---------------------------------------------------------------------------- #
agenaPats <- as.character(unique(agenaData$Individual.name))
columns <- c("ID", "Father", "Mother", "sex", "famid", "HLX", "KALRN", "CLSTN2", "CLSTN2_2", "LOC107986324", "SEMA6A", "TNKS", "TNKS_2", "NR6A1", "CALM1", "TENM1", "LRRK2", "CORO1C", "GRIK4", "RAB2B1M", "RAB2B1R", "RAB2B2", "Affected", "Status", "Time", "Subject", "LRRK2Clin", "OtherMutFL", "MutIsHomFL")
# 1, 2, 3, 4, 5, 6, 7, 8, 9, 
#LRRK2 = 12 in LUT, 
inputData = data.frame(matrix(nrow=0, ncol=length(columns)))
colnames(inputData) = columns
idProb = data.frame(matrix(nrow=0, ncol=1))

for (i in agenaPats) {
  subPersData <- subset(agenaData, agenaData$Individual.name == i)
  if (any(subPersData$Assay == assayLabels[12])) { # & any(subPersData$Assay == assayLabels[2])) {
    if (any(pedData$V9 == i)) {
      pedIdx <- which(pedData$V9 == i, arr.ind = TRUE)
      if (pedData$V3[pedIdx[1]] == 1 & pedData$V4[pedIdx[1]] == 1 & pedData$V6[pedIdx[1]] == "") {
        pedData$V3[pedIdx] <- 0
        pedData$V4[pedIdx] <- 0
        if (substr(pedData$V9[pedIdx[1]], 8, 8) == "_") {
          ssTemp <- matrix(unlist(strsplit(as.character(pedData$V9[pedIdx[1]]), split="_")), ncol=2)
          inputData[nrow(inputData)+1, 1] <- as.character(ssTemp[2])
          inputData[nrow(inputData), 5] <- as.character(ssTemp[1])
          inputData[nrow(inputData), 26] <- as.character(pedData$V9[pedIdx[1]])
        } else {
          inputData[nrow(inputData)+1, 1] <- as.character(pedData$V9[pedIdx], as.character("1"), sep="_")
          #inputData[nrow(inputData)+1, 1] <- as.character("1")
          inputData[nrow(inputData), 5] <- as.character(pedData$V9[pedIdx])
          inputData[nrow(inputData), 26] <- paste(as.character(pedData$V9[pedIdx]), as.character("1"), sep="_")
        }
      } else if (pedData$V1[pedIdx] == "") {
        if (substr(pedData$V9[pedIdx[1]], 8, 8) == "_") {
          ssTemp <- matrix(unlist(strsplit(as.character(pedData$V9[pedIdx[1]]), split="_")), ncol=2)
          inputData[nrow(inputData)+1, 1] <- as.character(ssTemp[2])
          inputData[nrow(inputData), 5] <- as.character(ssTemp[1])
          inputData[nrow(inputData), 26] <- as.character(pedData$V9[pedIdx[1]])
        } else {
          inputData[nrow(inputData)+1, 1] <- paste(as.character(pedData$V9[pedIdx]), as.character("1"), sep="_")
          #inputData[nrow(inputData)+1, 1] <- as.character("1")
          inputData[nrow(inputData), 5] <- as.character(pedData$V9[pedIdx])
          inputData[nrow(inputData), 26] <- paste(as.character(pedData$V9[pedIdx]), as.character("1"), sep="_")
        }
      } else {
        inputData[nrow(inputData)+1, 1] <- as.character(pedData$V2[pedIdx])
        inputData[nrow(inputData), 5] <- as.character(pedData$V1[pedIdx])
        inputData[nrow(inputData), 26] <- as.character(pedData$V9[pedIdx])
      }
      
      if (pedData$V3[pedIdx[1]] == 0) {
        inputData[nrow(inputData), 2] <- pedData$V3[pedIdx]
      } else {
        idxPop <- which(agenaPats == paste(inputData[nrow(inputData), 5], pedData$V3[pedIdx], sep="_"))
        if (length(idxPop) == 0) {
          inputData[nrow(inputData), 2] <- 0
        } else {
          inputData[nrow(inputData), 2] <- paste(pedData$V6[pedIdx], pedData$V3[pedIdx], sep="_")
        }
      }
      if (pedData$V4[pedIdx[1]] == 0) {
        inputData[nrow(inputData), 3] <- pedData$V4[pedIdx]
      } else {
        idxMom <- which(agenaPats == paste(inputData[nrow(inputData), 5], pedData$V4[pedIdx], sep="_"))
        if (length(idxMom) == 0) {
          inputData[nrow(inputData), 3] <- 0
        } else {
          inputData[nrow(inputData), 3] <- paste(pedData$V6[pedIdx], pedData$V4[pedIdx], sep="_")
        }
      }
      
      if (pedData$V10[pedIdx] == "F") {
        inputData[nrow(inputData),4] <- 2
      } else if (pedData$V10[pedIdx] == "M") {
        inputData[nrow(inputData),4] <- 1
      } else if (pedData$V10[pedIdx] == "U") {
        inputData[nrow(inputData),4] <- 3
      } else {
        print("Invalid gender detected")
      }
      
      sg <- data.frame(matrix(nrow=1, ncol=length(assayLabels)))
      for (k in 1:length(assayLabels)) {
        sg[k] <- subPersData$Replicates.Result[which(subPersData$Assay == assayLabels[k])]
        
        if (sg[k] == assayWT[k]) {
          inputData[nrow(inputData), 5+k] <- 0
          assayCounts[k, 1] <- assayCounts[k, 1] + 2
          assayCounts[k, 3] <- assayCounts[k, 3] + 1
          assayCounts[k, 6] <- assayCounts[k, 6] + 1
        } else if (sg[k] == assayHet[k]) {
          inputData[nrow(inputData), 5+k] <- 1
          assayCounts[k, 1] <- assayCounts[k, 1] + 1
          assayCounts[k, 2] <- assayCounts[k, 2] + 1
          assayCounts[k, 4] <- assayCounts[k, 4] + 1
          assayCounts[k, 6] <- assayCounts[k, 6] + 1
        } else if (sg[k] == assayHom[k]) {
          inputData[nrow(inputData), 5+k] <- 2
          assayCounts[k, 2] <- assayCounts[k, 2] + 2
          assayCounts[k, 5] <- assayCounts[k, 5] + 1
          assayCounts[k, 6] <- assayCounts[k, 6] + 1
        } else if (sg[k] == "mismatched" | sg[k] == "none") {
          inputData[nrow(inputData), 5+k] <- NaN
        } else {
          print(paste("ERROR: ", i, " - Uknown genotype for SNP-", as.character(k), " observed", sep = ""))
          if (sg[k] == "TC"){
            assayCounts[k, 1] <- assayCounts[k, 1] + 1
            assayCounts[k, 6] <- assayCounts[k, 6] + 1
          }
          inputData[nrow(inputData), 5+k] <- 0
        }
      }
      
      # Patient Diagnosis [affected] (0 = no diagnosis, 1 = PD)
      # Other tags to determine: Epilepsy; Depression
      aff <- subPersData$Primary.Current.Diagnosis[1]
      if (!is.na(aff)) {
        if (grepl("PD", aff)) {   # | grepl("ET", aff)) {
          inputData[nrow(inputData), 23] <- 1
        } else if (aff == "NA" | aff == "") {
          inputData[nrow(inputData), 23] <- 0
        } else {
          inputData[nrow(inputData), 23] <- 0
          print(paste("ERROR: Uknown affected state observed -", i, sep = " "))
        }
      } else {
        inputData[nrow(inputData), 23] <- 0
      }
      
      # Patient Status (0 = alive / censored, 1 = dead)
      inputData[nrow(inputData), 24] <- subPersData$Dead[1]
      
      # Patient TIME variable. If the patient was never diagnosed with PD, TIME is their current age
      # If the patient has PD, then TIME is age at first symptom presentation
      if (inputData[nrow(inputData), 23] == 1) {
        if (!is.na(subPersData$Age.at.Initial.Symptom[1])) {
          inputData[nrow(inputData), 25] <- subPersData$Age.at.Initial.Symptom[1]
        } else {
          
          inputData[nrow(inputData), 25] <- subPersData$Age.at.Initial.Symptom[1]
        }
      } else if (inputData[nrow(inputData), 23] == 0) {
        dobStr <- matrix(unlist(strsplit(as.character(subPersData$Date.of.Birth[1]), split="/")), ncol=3)
        if (!is.na(dobStr[1])) {
          if (strtoi(dobStr[3]) > 2008) {
            inputData[nrow(inputData), 25] <- 2008 - (strtoi(dobStr[3])-100)
          } else {
            inputData[nrow(inputData), 25] <- 2008 - strtoi(dobStr[3])
          }
        } else {
          inputData[nrow(inputData), 25] <- NaN
        }
      }
      
      if (subPersData$Mutation[1] != "") { #Mutation vector is not empty
        if (subPersData$Mutation[1] == lrrk2HetTag) {
          inputData[nrow(inputData), 27] <- 1
          inputData[nrow(inputData), 28] <- 0
          inputData[nrow(inputData), 29] <- 0
        } else if (subPersData$Mutation[1] == lrrk2HomTag) {
          inputData[nrow(inputData), 27] <- 1
          inputData[nrow(inputData), 28] <- 0
          inputData[nrow(inputData), 29] <- 0
        } else if (grepl(subPersData$Mutation[1], lrrk2HetTag)) {
          inputData[nrow(inputData), 27] <- 1
          inputData[nrow(inputData), 28] <- 1
          mutStr <- str_remove(subPersData$Mutation[1], lrrk2HetTag)
          if (grepl(mutString, "+/+") | grepl(mutString, "PINK1 Deletion")) {
            inputData[nrow(inputData), 29] <- 1
          } else {
            inputData[nrow(inputData), 29] <- 0
          }
        } else if (grepl(subPersData$Mutation[1], lrrk2HomTag)) {
          inputData[nrow(inputData), 27] <- 1
          inputData[nrow(inputData), 28] <- 1
          mutStr <- str_remove(subPersData$Mutation[1], lrrk2HetTag)
          if (grepl(mutString, "+/+") | grepl(mutString, "PINK1 Deletion")) {
            inputData[nrow(inputData), 29] <- 1
          } else {
            inputData[nrow(inputData), 29] <- 0
          }
        } else {
          inputData[nrow(inputData), 27] <- 0
          inputData[nrow(inputData), 28] <- 1
          if (grepl(subPersData$Mutation[1], "+/+") | grepl(subPersData$Mutation[1], "PINK1 Deletion")) {
            inputData[nrow(inputData), 29] <- 1
          } else {
            inputData[nrow(inputData), 29] <- 0
          }
        }
      } else { #Mutation vector is empty
        inputData[nrow(inputData), 27] <- 0
        inputData[nrow(inputData), 28] <- 0
        inputData[nrow(inputData), 29] <- 0
      }
      
      
      if (inputData[nrow(inputData),2] == 0 & inputData[nrow(inputData),3] != 0) {
        idProb[nrow(idProb)+1, ] = as.character(
          paste(as.character(inputData[nrow(inputData),5]), 
                as.character(inputData[nrow(inputData),3]),
                sep="_"))
      } else if (inputData[nrow(inputData),2] != 0 & inputData[nrow(inputData),3] == 0) {
        idProb[nrow(idProb)+1, ] = as.character(
          paste(as.character(inputData[nrow(inputData),5]), 
                as.character(inputData[nrow(inputData),2]),
                sep="_"))
      }
        
    }
  } else {
    idProb[nrow(idProb)+1, ] <- as.character(i)
  }
}

# ---------------------------------------------------------------------- #
# THIS SECTION IS WHERE WE SCAN THROUGH ALL FAILED SAMPLES TO MAKE SURE  #
# THAT THEY ARE ALSO REMOVED FROM THE MOM AND DAD ID COLUMNS             #
# ANYWHERE A REMOVED ID IS DETECTED, THE VALUE IS REPLACED WITH A 0      #
# ---------------------------------------------------------------------- #
# THIS SECTION HAS BEEN REDACTED
# Previously, we used this section to remove potentially problematic IDs
# from the inputData data frame because we would use that same data frame 
# to generate our kinship matrix. However, we now build the kinship
# matrix independently from our inputData data frame
#for (i in idProb[,1]) {
#  ssTemp <- matrix(unlist(strsplit(as.character(i), split="_")), ncol=2)
#  idxMomBad <- which(inputData$Mother == ssTemp[2] & inputData$famid == ssTemp[1])
#  if (length(idxMomBad) > 0) {
#    inputData$Mother[idxMomBad] <- 0
#    inputData$Father[idxMomBad] <- 0
#  }
#  idxDadBad <- which(inputData$Father == ssTemp[2] & inputData$famid == ssTemp[1])
#  if (length(idxDadBad) > 0) {
#    inputData$Father[idxDadBad] <- 0
#    inputData$Mother[idxDadBad] <- 0
#  }
#}

#ssttM <- setdiff(inputData$Mother, inputData$Subject)
#for (i in 2:length(ssttM)) {
#  idxDel <- which(inputData$Mother == ssttM[i])
#  inputData$Mother[idxDel] <- 0
#  inputData$Father[idxDel] <- 0
#}

#ssttF <- setdiff(inputData$Father, inputData$Subject)
#for (i in 2:length(ssttF)) {
#  idxDel <- which(inputData$Father == ssttF[i])
#  inputData$Father[idxDel] <- 0
#  inputData$Mother[idxDel] <- 0
#}

#idxDel2 <- which(inputData$Father == 0)
#inputData$Mother[idxDel2] <- 0

#idxDel2 <- which(inputData$Mother == 0)
#inputData$Father[idxDel2] <- 0
# ---------------------------------------------------------------------------- #
# THIS SECTION IS WHERE WE HARD CODE FIXES FOR KNOWN SIGNIFICANT MIXUP SAMPLES #
# ---------------------------------------------------------------------------- #

kinPats <- as.character(unique(pedData$V9))
columnsKin <- c("ID", "Father", "Mother", "sex", "famid")
inputKin = data.frame(matrix(nrow=0, ncol=length(columnsKin)))
colnames(inputKin) = columnsKin
for (i in kinPats) {
  pedIdx <- which(pedData$V9 == i, arr.ind = TRUE)
  if (pedData$V3[pedIdx] == 1 & pedData$V4[pedIdx] == 1 & pedData$V6[pedIdx] == "") {
    pedData$V3[pedIdx] <- 0
    pedData$V4[pedIdx] <- 0
    
    if (substr(pedData$V9[pedIdx], 8, 8) == "_") {
      ssTemp <- matrix(unlist(strsplit(as.character(pedData$V9[pedIdx]), split="_")), ncol=2)
      inputKin[nrow(inputKin)+1, 1] <- as.character(pedData$V9[pedIdx])
      inputKin[nrow(inputKin), 5] <- as.character(ssTemp[1])
    } else {
      inputKin[nrow(inputKin)+1, 1] <- paste(as.character(pedData$V9[pedIdx]), as.character("1"), sep="_")
      inputKin[nrow(inputKin), 5] <- as.character(pedData$V9[pedIdx])
    }
  } else if (pedData$V1[pedIdx] == "") {
    if (substr(pedData$V9[pedIdx], 8, 8) == "_") {
      ssTemp <- matrix(unlist(strsplit(as.character(pedData$V9[pedIdx[1]]), split="_")), ncol=2)
      inputKin[nrow(inputKin)+1, 1] <- as.character(pedData$V9[pedIdx[1]])
      inputKin[nrow(inputKin), 5] <- as.character(ssTemp[1])
    } else {
      inputKin[nrow(inputKin)+1, 1] <- paste(as.character(pedData$V9[pedIdx]), as.character("1"), sep="_")
      #inputData[nrow(inputData)+1, 1] <- as.character("1")
      inputKin[nrow(inputKin), 5] <- as.character(pedData$V9[pedIdx])
    }
  } else {
    inputKin[nrow(inputKin)+1, 1] <- as.character(pedData$V9[pedIdx])
    inputKin[nrow(inputKin), 5] <- as.character(pedData$V1[pedIdx])
  }
  
  if (pedData$V3[pedIdx] == 0) {
    inputKin[nrow(inputKin), 2] <- pedData$V3[pedIdx]
  } else {
    idxPop <- which(kinPats == paste(inputKin[nrow(inputKin), 5], pedData$V3[pedIdx], sep="_"))
    if (length(idxPop) == 0) {
      inputKin[nrow(inputKin), 2] <- 0
    } else {
      inputKin[nrow(inputKin), 2] <- paste(pedData$V6[pedIdx], pedData$V3[pedIdx], sep="_")
    }
  }
  if (pedData$V4[pedIdx] == 0) {
    inputKin[nrow(inputKin), 3] <- pedData$V4[pedIdx]
  } else {
    idxMom <- which(agenaPats == paste(inputKin[nrow(inputKin), 5], pedData$V4[pedIdx], sep="_"))
    if (length(idxMom) == 0) {
      inputKin[nrow(inputKin), 3] <- 0
    } else {
      inputKin[nrow(inputKin), 3] <- paste(pedData$V6[pedIdx], pedData$V4[pedIdx], sep="_")
    }
  }
  
  if (pedData$V10[pedIdx] == "F") {
    inputKin[nrow(inputKin),4] <- 2
  } else if (pedData$V10[pedIdx] == "M") {
    inputKin[nrow(inputKin),4] <- 1
  } else if (pedData$V10[pedIdx] == "U") {
    inputKin[nrow(inputKin),4] <- 3
  } else {
    print("Invalid gender detected")
  }
}

idxNoDadKin <- which(inputKin$Father == 0 & inputKin$Mother != 0, arr.ind = TRUE)
inputKin$Mother[idxNoDadKin] <- 0

idxNoMomKin <- which(inputKin$Mother == 0 & inputKin$Father != 0, arr.ind = TRUE)
inputKin$Father[idxNoMomKin] <- 0

inputData_Affected <- subset(inputData, inputData$Affected == 1)
inputData_AffAndG2019S <- subset(inputData_Affected, inputData_Affected$LRRK2 == 1)

pedTree <- pedigree(inputKin$ID, inputKin$Father, inputKin$Mother, sex=inputKin$sex, famid=inputKin$famid, missid="0")
kinshipMatrix <- kinship(pedTree)
#print(pedTree[3])
#plot(pedTree[3])

inputData_NoOtherMuts <- subset(inputData, inputData$OtherMutFL == 0) # & !(inputData$LRRK2 == 0 & inputData$OtherMutFL == 1 & inputData$Affected == 1))
inputData_NoHomozMuts <- subset(inputData, !(inputData$LRRK2 == 0 & inputData$MutIsHomFL == 1 & inputData$Affected == 1))
inputData_NoHomozMuts <- subset(inputData, !(inputData$LRRK2 == 0 & inputData$OtherMutFL == 1 & inputData$Affected == 1))
rejectedData_OtherMut <- subset(inputData, inputData$OtherMutFL == 1)
rejectedData_OtherMut$MutationNotes <- c()
for (i in 1:length(rejectedData_OtherMut$Subject)) {
  nameOI <- rejectedData_OtherMut$Subject[i]
  if (substr(nameOI, start=5, stop=8) == "case") {
    nameOI = substr(nameOI, start=1, stop=nchar(nameOI)-2)
  }
  subPersData <- subset(agenaData, agenaData$Individual.name == nameOI)
  rejectedData_OtherMut$MutationNotes[i] <- subPersData$Mutation[1]
}
write_xlsx(inputData_NoOtherMuts, "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\inputData_NoMutationStrictCohort_Updated-08-25-23.xlsx")
write_xlsx(inputData_NoOtherMuts, "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\inputData_NoMutationLiberalCohort_Updated-08-25-23.xlsx")

write_xlsx(inputData_No)
write_xlsx(rejectedData_OtherMut, "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\OmittedPatients_NoMutationCohort_Updated-08-25-23.xlsx")

strictNoMutCohortPath <- "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\inputData_NoMutationStrictCohort_Updated-08-25-23.xlsx"
liberalNoMutCohortPath <- "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\inputData_NoMutationLiberalCohort_Updated-08-25-23.xlsx"

inputData_NoOtherMuts <- read_excel(liberalNoMutCohortPath)

# INPUT DATASET METRICS (COHORT)
#  First lets do MAF so that we can change our homozygous 2 calls
# to 1's, then change 2s to 1s at the end.

geneCols <- c("HLX", "KALRN", "CLSTN2", "CLSTN2_2", "LOC107986324", "SEMA6A", "TNKS", "TNKS_2", "NR6A1", "CALM1", "TENM1", "LRRK2", "CORO1C", "GRIK4", "RAB2B1M", "RAB2B1R", "RAB2B2")
for (gc in geneCols) {
  print(gc)
  print(length(which(!is.na(inputData_NoOtherMuts[[gc]]))))
  print((length(which(inputData_NoOtherMuts[[gc]] == 0)) * 2) + (length(which(inputData_NoOtherMuts[[gc]] == 1)) * 1))
  print((length(which(inputData_NoOtherMuts[[gc]] == 2)) * 2) + (length(which(inputData_NoOtherMuts[[gc]] == 1)) * 1))
  print(length(which(inputData_NoOtherMuts[[gc]] == 0)))
  print(length(which(inputData_NoOtherMuts[[gc]] == 1)))
  print(length(which(inputData_NoOtherMuts[[gc]] == 2)))
  
  inputData_NoOtherMuts[[gc]][which(inputData_NoOtherMuts[[gc]] == 2)] <- 1
}

#Number of families
print()
#Number of Males
print(length(which(inputData_NoOtherMuts$sex == 1)))
#Number of Females
print(length(which(inputData_NoOtherMuts$sex == 2)))
#Number undeclared gender
print(length(which(inputData_NoOtherMuts$sex == 3)))
#Number GS+ Carriers
print(length(which(inputData_NoOtherMuts$LRRK2 == 1)))
#Number affected PD
print(length(which(inputData_NoOtherMuts$Affected == 1)))
#Number affected PD & GS+
print(length(which(inputData_NoOtherMuts$Affected == 1 & inputData_NoOtherMuts$LRRK2 == 1)))
#Average AAO
print(mean(inputData_NoOtherMuts$Time[!is.na(inputData_NoOtherMuts$Time)]))

# "HLX", "KALRN", "CLSTN2", "CLSTN2_2", "LOC107986324", "SEMA6A", "TNKS", "TNKS_2", "NR6A1", 
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + HLX + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + KALRN + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + CLSTN2 + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + CLSTN2_2 + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + LOC107986324 + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + SEMA6A + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + TNKS + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + TNKS_2 + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + NR6A1 + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
# "CALM1", "TENM1", "LRRK2", "CORO1C", "GRIK4", "RAB2B1M", "RAB2B1R", "RAB2B2"
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + CALM1 + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + TENM1 + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + CORO1C + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + GRIK4 + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + RAB2B1M + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + RAB2B1R + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)
penetranceModel_NoOtherMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + RAB2B2 + (1|Subject), data=inputData_NoOtherMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoOtherMuts)
confint(penetranceModel_NoOtherMuts, level=0.95)

inputDataAff_NoOtherMuts <- subset(inputData_NoOtherMuts, inputData_NoOtherMuts$Affected == 1)
for (gc in geneCols) {
  # index of all affected individuals with LRRK2+ and SNPoi+ 
  idxOIpp <- which(inputDataAff_NoOtherMuts$LRRK2==1 & inputDataAff_NoOtherMuts[[gc]]==1 & !is.na(inputDataAff_NoOtherMuts$Time))
  # index of all affected individuals with LRRK2+ and SNPoi-
  idxOIpn <- which(inputDataAff_NoOtherMuts$LRRK2==1 & inputDataAff_NoOtherMuts[[gc]]==0 & !is.na(inputDataAff_NoOtherMuts$Time))
  
  # 1- SNPoi Name
  print(gc)
  # 5- Mean AAO for LRRK2+ and SNPoi+
  print(mean(inputDataAff_NoOtherMuts$Time[idxOIpp]))
  # 6- Standard deviation AAO for LRRK+ and SNPoi+
  print(sd(inputDataAff_NoOtherMuts$Time[idxOIpp]))
  # 7- N for LRRK2+ and SNPoi+
  print(length(idxOIpp))
  # 8- Mean AAO for LRRK2+ and SNPoi+
  print(mean(inputDataAff_NoOtherMuts$Time[idxOIpn]))
  # 9- Standard deviation AAO for LRRK+ and SNPoi+
  print(sd(inputDataAff_NoOtherMuts$Time[idxOIpn]))
  # 10- N for LRRK2+ and SNPoi+
  print(length(idxOIpn))
  
  # 2- Number of SNPoi carriers
  print(length(which(inputDataAff_NoOtherMuts[[gc]] == 1 & !is.na(inputDataAff_NoOtherMuts$Time))))
  # 3- Number SNPoi data missing
  print(length(which(is.na(inputDataAff_NoOtherMuts[[gc]]))))
  # 4- Number of SNPoi carriers with missing time
  print(length(which(inputDataAff_NoOtherMuts[[gc]] == 1 & is.na(inputDataAff_NoOtherMuts$Time))))
}

# "HLX", "KALRN", "CLSTN2", "CLSTN2_2", "LOC107986324", "SEMA6A", "TNKS", "TNKS_2", "NR6A1", 
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + HLX + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + KALRN + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + CLSTN2 + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + CLSTN2_2 + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + LOC107986324 + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + SEMA6A + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + TNKS + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + TNKS_2 + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + NR6A1 + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
# "CALM1", "TENM1", "LRRK2", "CORO1C", "GRIK4", "RAB2B1M", "RAB2B1R", "RAB2B2"
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + CALM1 + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + TENM1 + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + CORO1C + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + GRIK4 + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + RAB2B1M + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + RAB2B1R + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)
ageOnsetModel_NoOtherMuts <- lmekin(Time ~ sex + LRRK2 + RAB2B2 + (1|Subject), data=inputDataAff_NoOtherMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoOtherMuts)


# ------------ No homozygous mutations ---------------- #
for (gc in geneCols) {
  print(gc)
  print(length(which(!is.na(inputData_NoHomozMuts[[gc]]))))
  print((length(which(inputData_NoHomozMuts[[gc]] == 0)) * 2) + (length(which(inputData_NoHomozMuts[[gc]] == 1)) * 1))
  print((length(which(inputData_NoHomozMuts[[gc]] == 2)) * 2) + (length(which(inputData_NoHomozMuts[[gc]] == 1)) * 1))
  print(length(which(inputData_NoHomozMuts[[gc]] == 0)))
  print(length(which(inputData_NoHomozMuts[[gc]] == 1)))
  print(length(which(inputData_NoHomozMuts[[gc]] == 2)))
  
  inputData_NoHomozMuts[[gc]][which(inputData_NoHomozMuts[[gc]] == 2)] <- 1
}

#Number of families
print()
#Number of Males
print(length(which(inputData_NoHomozMuts$sex == 1)))
#Number of Females
print(length(which(inputData_NoHomozMuts$sex == 2)))
#Number undeclared gender
print(length(which(inputData_NoHomozMuts$sex == 3)))
#Number GS+ Carriers
print(length(which(inputData_NoHomozMuts$LRRK2 == 1)))
#Number affected PD
print(length(which(inputData_NoHomozMuts$Affected == 1)))
#Number affected PD & GS+
print(length(which(inputData_NoHomozMuts$Affected == 1 & inputData_NoHomozMuts$LRRK2 == 1)))
#Average AAO
print(mean(inputData_NoHomozMuts$Time[!is.na(inputData_NoHomozMuts$Time)]))

# "HLX", "KALRN", "CLSTN2", "CLSTN2_2", "LOC107986324", "SEMA6A", "TNKS", "TNKS_2", "NR6A1", 
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + HLX + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + KALRN + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + CLSTN2 + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + CLSTN2_2 + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + LOC107986324 + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + SEMA6A + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + TNKS + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + TNKS_2 + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + NR6A1 + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
# "CALM1", "TENM1", "LRRK2", "CORO1C", "GRIK4", "RAB2B1M", "RAB2B1R", "RAB2B2"
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + CALM1 + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + TENM1 + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + CORO1C + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + GRIK4 + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + RAB2B1M + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + RAB2B1R + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)
penetranceModel_NoHomozMuts <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + RAB2B2 + (1|Subject), data=inputData_NoHomozMuts, varlist=2*kinshipMatrix)
print(penetranceModel_NoHomozMuts)

inputDataAff_NoHomozMuts <- subset(inputData_NoHomozMuts, inputData_NoHomozMuts$Affected == 1)
for (gc in geneCols) {
  # index of all affected individuals with LRRK2+ and SNPoi+ 
  idxOIpp <- which(inputDataAff_NoHomozMuts$LRRK2==1 & inputDataAff_NoHomozMuts[[gc]]==1 & !is.na(inputDataAff_NoHomozMuts$Time))
  # index of all affected individuals with LRRK2+ and SNPoi-
  idxOIpn <- which(inputDataAff_NoHomozMuts$LRRK2==1 & inputDataAff_NoHomozMuts[[gc]]==0 & !is.na(inputDataAff_NoHomozMuts$Time))
  
  # 1- SNPoi Name
  print(gc)
  # 5- Mean AAO for LRRK2+ and SNPoi+
  print(mean(inputDataAff_NoHomozMuts$Time[idxOIpp]))
  # 6- Standard deviation AAO for LRRK+ and SNPoi+
  print(sd(inputDataAff_NoHomozMuts$Time[idxOIpp]))
  # 7- N for LRRK2+ and SNPoi+
  print(length(idxOIpp))
  # 8- Mean AAO for LRRK2+ and SNPoi+
  print(mean(inputDataAff_NoHomozMuts$Time[idxOIpn]))
  # 9- Standard deviation AAO for LRRK+ and SNPoi+
  print(sd(inputDataAff_NoHomozMuts$Time[idxOIpn]))
  # 10- N for LRRK2+ and SNPoi+
  print(length(idxOIpn))
  
  # 2- Number of SNPoi carriers
  print(length(which(inputDataAff_NoHomozMuts[[gc]] == 1 & !is.na(inputDataAff_NoHomozMuts$Time))))
  # 3- Number SNPoi data missing
  print(length(which(is.na(inputDataAff_NoHomozMuts[[gc]]))))
  # 4- Number of SNPoi carriers with missing time
  print(length(which(inputDataAff_NoHomozMuts[[gc]] == 1 & is.na(inputDataAff_NoHomozMuts$Time))))
  
  #print(length(which(inputDataAff_NoHomozMuts[[gc]] == 1 & !is.na(inputDataAff_NoHomozMuts$Time))))
}

# "HLX", "KALRN", "CLSTN2", "CLSTN2_2", "LOC107986324", "SEMA6A", "TNKS", "TNKS_2", "NR6A1", 
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + HLX + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + KALRN + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + CLSTN2 + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + CLSTN2_2 + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + LOC107986324 + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + SEMA6A + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + TNKS + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + TNKS_2 + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + NR6A1 + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
# "CALM1", "TENM1", "LRRK2", "CORO1C", "GRIK4", "RAB2B1M", "RAB2B1R", "RAB2B2"
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + CALM1 + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + TENM1 + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + CORO1C + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + GRIK4 + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + RAB2B1M + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + RAB2B1R + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)
ageOnsetModel_NoHomozMuts <- lmekin(Time ~ sex + LRRK2 + RAB2B2 + (1|Subject), data=inputDataAff_NoHomozMuts, varlist=2*kinshipMatrix)
print(ageOnsetModel_NoHomozMuts)

penetranceModelAll <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + HLX + (1|Subject), data=inputData, varlist=2*kinshipMatrix)
penetranceModelAff <- coxme(Surv(Time, Affected) ~ sex + LRRK2 + CORO1C + (1|Subject), data=inputData_Affected, varlist=2*kinshipMatrix)
penetranceModelG2S <- coxme(Surv(Time, Affected) ~ sex + CORO1C + (1|Subject), data=inputData_G2019S, varlist=2*kinshipMatrix)
print(penetranceModelAll)
print(penetranceModelAff)
print(penetranceModelG2S)

penetranceModelPH <- coxph(Surv(Time, Affected) ~ sex + LRRK2 + HLX, data=inputData)
print(penetranceModelPH)


ageOnsetModel <- lmekin(Time ~ sex + LRRK2 + HLX + (1|Subject), data=inputData_Affected, varlist=2*kinshipMatrix)
print(ageOnsetModel)


# SUBSET ONLY AFFECTED INDIVIDUALS #
inputData_CORO1C <- subset(inputData, !is.nan(inputData$CORO1C)) #inputData$CORO1C!=2 & 
inputData_CORO1C$CORO1C[inputData_CORO1C$CORO1C == 0] = "Wildtype"
inputData_CORO1C$CORO1C[inputData_CORO1C$CORO1C == 1] = "Mutant"

#compriskEst <- tidycmprsk::crr(Surv(Time, as.factor(Affected)) ~ CORO1C, data=inputData_CORO1C)
compriskFit <- tidycmprsk::cuminc(Surv(Time, as.factor(Affected)) ~ CORO1C, data=inputData_CORO1C)
print(compriskFit)

compriskFit %>%
  ggcuminc() + 
  add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})",
                stats_label = list(n.risk = "Number at Risk", cum.event = "Number Affected")) + 
  scale_color_manual(name="Genotype", values=c("red", "blue"), guide="legend") +
  scale_fill_manual(name="Genotype", values=c("red", "blue"), guide="legend") +
  labs(title="Coro1c on Risk", x="Time (Years)", y="Cumulative Incidence") + 
  theme(legend.position="bottom") +
  theme_minimal()

# SUBSET ONLY AFFECTED INDIVIDUALS #
inputData_RAB2B1M <- subset(inputData, !is.nan(inputData$RAB2B1M)) #inputData$RAB2B1M!=2 & 
inputData_RAB2B1M$RAB2B1M[inputData_RAB2B1M$RAB2B1M == 0] = "Wildtype"
inputData_RAB2B1M$RAB2B1M[inputData_RAB2B1M$RAB2B1M == 1] = "Mutant"

#compriskEst <- tidycmprsk::crr(Surv(Time, as.factor(Affected)) ~ CORO1C, data=inputData_CORO1C)
compriskFit <- tidycmprsk::cuminc(Surv(Time, as.factor(Affected)) ~ RAB2B1M, data=inputData_RAB2B1M)
print(compriskFit)

compriskFit %>%
  ggcuminc() + 
  add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})",
                stats_label = list(n.risk = "Number at Risk", cum.event = "Number Affected")) + 
  scale_color_manual(name="Genotype", values=c("red", "blue"), guide="legend") +
  scale_fill_manual(name="Genotype", values=c("red", "blue"), guide="legend") +
  labs(title="RAB2B1M on Risk", x="Time (Years)", y="Cumulative Incidence") + 
  theme(legend.position="bottom") +
  theme_minimal()


inputData_RAB2B1R <- subset(inputData, !is.nan(inputData$RAB2B1R)) #inputData$RAB2B1M!=2 & 
inputData_RAB2B1R$RAB2B1R[inputData_RAB2B1R$RAB2B1R == 0] = "Wildtype"
inputData_RAB2B1R$RAB2B1R[inputData_RAB2B1R$RAB2B1R == 1] = "Mutant"

#compriskEst <- tidycmprsk::crr(Surv(Time, as.factor(Affected)) ~ CORO1C, data=inputData_CORO1C)
compriskFit <- tidycmprsk::cuminc(Surv(Time, as.factor(Affected)) ~ RAB2B1R, data=inputData_RAB2B1R)
print(compriskFit)

compriskFit %>%
  ggcuminc() + 
  add_confidence_interval() +
  add_risktable(risktable_stats = "{n.risk} ({cum.event})",
                stats_label = list(n.risk = "Number at Risk", cum.event = "Number Affected")) + 
  scale_color_manual(name="Genotype", values=c("red", "blue"), guide="legend") +
  scale_fill_manual(name="Genotype", values=c("red", "blue"), guide="legend") +
  labs(title="RAB2B1R on Risk", x="Time (Years)", y="Cumulative Incidence") + 
  theme(legend.position="bottom") +
  theme_minimal()

ccLP_GoodOnly <-  Cuminc(time="Time", status="Affected", group="Geno2", data=inputData_GoodOnly)
ccLP_GoodOnlyP2 <- cuminc(ftime=inputData_GoodOnly$Time, fstatus=inputData_GoodOnly$Affected, group=inputData_GoodOnly$Geno2)
plot(ccLP_GoodOnlyP2, xlab="Years", ylab="Probability", xlim=c(0,100), ylim=(c(0,1)))



inputData_GoodOnly2 <- subset(inputData, inputData$Geno2!=2 & !is.nan(inputData$Geno2))
inputData_AffOnly <- subset(inputData, inputData$Affected==1 & !is.nan(inputData$Time) & !is.nan(inputData$Geno2))
ccLP_GoodOnly2 <-  cuminc(ftime=inputData_GoodOnly2$Time, fstatus=inputData_GoodOnly2$Affected, group=inputData_GoodOnly2$Geno2)
plot(ccLP_GoodOnly2, xlab="Years", ylab="Probability", xlim=c(0,100), ylim=(c(0,1)))


inputData_GSOnly <- subset(inputData, (inputData$Geno1==1 | inputData$Geno1==2) & !is.nan(inputData$Geno2))
ccLP_gsOnly <- cuminc(ftime=inputData_GSOnly$Time, fstatus=inputData_GSOnly$Affected, group=inputData_GSOnly$Geno2)
plot(ccLP_gsOnly, xlab="Years", ylab="Probability", xlim=c(0,100), ylim=(c(0,1)))

ccLP <- cuminc(Surv(inputData$Time, inputData$Affected) ~ inputData$Geno2)
ccLP <- cuminc(ftime=inputData$Time, fstatus=inputData$Affected, group=inputData$Geno2)
plot(ccLP, xlab="Years", ylab="Probability", xlim=c(0,100), ylim=(c(0,1)))
ccLP2 <- cuminc(ftime=inputData_AffOnly$Time, fstatus=inputData_AffOnly$Affected, group=inputData_AffOnly$Geno2)
plot(ccLP2, xlab="Years", ylab="Probability", xlim=c(0,100), ylim=(c(0,1)))

ageonsetModel <- lmekin(Surv(Time, Affected) ~ sex + Geno1 + Geno2, data=inputData, varlist=2*kinshipMatrix)
