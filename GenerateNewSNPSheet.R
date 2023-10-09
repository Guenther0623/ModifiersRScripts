# This function was used to generate new Genotyping files.
# -First finds all unique IDs from the Agena SNP output
# -Fetch patient/clinical data from PED file by matching IDs
# -Match unique/Progeny IDs and data to Agena SNPs
# -Recompile, organize, and save

clinicalDIR <- "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\"
clinicalCSV <- "Info_Tunisia_samples.csv"
clinicalFUL <- paste(clinicalDIR, clinicalCSV, sep="")

agenasnpDIR <- "C:\\Users\\HAL2\\Desktop\\Tunisia Modifiers Data\\"
agenasnpCSV <- "08.08.23_Tunisia_All_Info_SNP_Calls_DylanEdit.csv"
agenasnpFUL <- paste(agenasnpDIR, agenasnpCSV, sep="")

clinicalData <- read.csv(clinicalFUL)
agenasnpData <- read.csv(agenasnpFUL)

geneIDs <- c("HLX,DUSP10", "KALRN", "KALRN", "CLSTN2", "CLSTN2", "RAP2B,ARHGEF26", "RAP2B_A", 
             "RAP2B,ARHGEF26", "RAP2B_B", "LOC107986324", "LOC107986324", "SEMA6A", "SEMA6A", "TNKS", "TNKS", 
             "NR6A1", "NR6A1", "GRIK4", "GRIK4", "LRRK2", "LRRK2", "CORO1C", "CORO1C", 
             "CALM1,TTC7B", "CALM1", "TENM1", "TENM1")
assayIDs <- c("rs141686162", "rs145611031", "rs114392312", "rs150382576", "rs140825591", "rs59679443", "rs73007728", 
              "rs16846845", "rs58923546", "rs12272007", "rs6824505", "rs73781088", "rs17509902", "rs28398294", "rs28511573", 
              "rs148922482", "rs150792641", "rs28470321", "rs10502241", "rs34637584", "rs190882284", "rs77395454", "rs3825252", 
              "rs76788674", "rs76039170", "rs185981774", "rs185078407")

MainReduLabels <- c("Main", "Main", "Redundant", "Main", "Redundant", "Main", "Redundant", 
                    "Main", "Redundant", "Main", "Redundant", "Main", "Redundant", "Main", "Redundant", 
                    "Main", "Redundant", "Main", "Redundant", "Main", "Redundant", "Main", "Redundant",
                    "Main", "Redundant", "Main", "Redundant")



agenaPats <- as.character(unique(agenasnpData$Individual.name))

outputFrameCols <- c("", "Assay", "Individual name", "Replicates Result", "Pedigree name", 
                     "Gender", "Date of Birth", "Dead", "Primary/Current Diagnosis", 
                     "Age at Initial Symptom", "Age at Primary Diagnosis", "Individual notes", 
                     "Mutation", "Gene", "Main or Redundant", "Was Manually Called?")

outputData = data.frame(matrix(nrow=0, ncol=length(outputFrameCols)))
colnames(outputData) = outputFrameCols
rowNum = 1

for (i in agenaPats) {
  clinicPatIDX <- which(clinicalData$Individual.name == i)
  if (!(length(clinicPatIDX) == 0)) {
    agenaPatIDX <- which(agenasnpData$Individual.name == i)
    
    possPedNames <- as.character(unique(clinicalData$Pedigree.name[clinicPatIDX]))
    possGenders <- as.character(unique(clinicalData$Gender[clinicPatIDX]))
    possBirthDates <- as.character(unique(clinicalData$Date.of.Birth[clinicPatIDX]))
    possDeadStates <- unique(clinicalData$Dead[clinicPatIDX])
    possPrimaryDiag <- as.character(unique(clinicalData$Primary.Current.Diagnosis[clinicPatIDX]))
    possAgeInitSym <- as.character(unique(clinicalData$Age.at.Initial.Symptom[clinicPatIDX]))
    possAgePrimDia <- as.character(unique(clinicalData$Age.at.Primary.Diagnosis[clinicPatIDX]))
    possIndivNotes <- as.character(unique(clinicalData$Individual.Notes[clinicPatIDX]))
    possMutations <- as.character(unique(clinicalData$Mutation[clinicPatIDX]))
    
    if (length(possPedNames) > 1) {
      print(paste("ERROR: Duplicate possibilities found for", i, sep=" "))
    }
    if (length(possGenders) > 1) {
      print(paste("ERROR: Duplicate possibilities found for", i, sep=" "))
    }
    if (length(possBirthDates) > 1) {
      print(paste("ERROR: Duplicate possibilities found for", i, sep=" "))
    }
    if (length(possDeadStates) > 1) {
      print(paste("ERROR: Duplicate possibilities found for", i, sep=" "))
    }
    if (length(possPrimaryDiag) > 1) {
      print(paste("ERROR: Duplicate possibilities found for", i, sep=" "))
    }
    if (length(possAgeInitSym) > 1) {
      print(paste("ERROR: Duplicate possibilities found for", i, sep=" "))
    }
    if (length(possAgePrimDia) > 1) {
      print(paste("ERROR: Duplicate possibilities found for", i, sep=" "))
    }
    if (length(possIndivNotes) > 1) {
      print(paste("ERROR: Duplicate possibilities found for", i, sep=" "))
    }
    if (length(possMutations) > 1) {
      print(paste("ERROR: Duplicate possibilities found for", i, sep=" "))
    }
    
    if (is.na(possPedNames)) {
      possPedNames <- ""
    }
    if (is.na(possGenders)) {
      possGenders <- ""
    }
    if (is.na(possBirthDates)) {
      possBirthDates <- ""
    }
    if (is.na(possDeadStates)) {
      possDeadStates <- ""
    }
    if (is.na(possPrimaryDiag)) {
      possPrimaryDiag <- ""
    }
    if (is.na(possAgeInitSym)) {
      possAgeInitSym <- ""
    }
    if (is.na(possAgePrimDia)) {
      possAgePrimDia <- ""
    }
    if (is.na(possIndivNotes)) {
      possIndivNotes <- ""
    }
    if (is.na(possMutations)) {
      possMutations <- ""
    }
    indivName <- i
    
    if (possBirthDates == "") {
      possBirthString = ""
    } else {
      birthSubStr <- strsplit(possBirthDates, split="/")
      possBirthString <- paste(birthSubStr[[1]][1], birthSubStr[[1]][2], paste("19", birthSubStr[[1]][3], sep=""), sep="/")
    }
    
    for (j in 1:length(geneIDs)) {
      iidx <- which(agenasnpData$Assay[agenaPatIDX] == assayIDs[j])
      outputData[rowNum, 1] <- rowNum
      outputData[rowNum, 2] <- assayIDs[j]
      outputData[rowNum, 3] <- indivName
      # 04
      outputData[rowNum, 5] <- possPedNames
      outputData[rowNum, 6] <- possGenders
      outputData[rowNum, 7] <- possBirthString
      outputData[rowNum, 8] <- possDeadStates
      outputData[rowNum, 9] <- possPrimaryDiag
      outputData[rowNum, 10] <- possAgeInitSym
      outputData[rowNum, 11] <- possAgePrimDia
      outputData[rowNum, 12] <- possIndivNotes
      outputData[rowNum, 13] <- possMutations
      outputData[rowNum, 14] <- geneIDs[j]
      outputData[rowNum, 15] <- MainReduLabels[j]
      # 16
      
      if (length(iidx) == 0) {
        outputData[rowNum, 4] <- "none"
        outputData[rowNum, 16] <- ""
      } else if (length(iidx) == 1) {
        outputData[rowNum, 4] <- agenasnpData$Replicates.Result[agenaPatIDX[iidx]]
        outputData[rowNum, 16] <- agenasnpData$Was.Manually.Called.[agenaPatIDX[iidx]]
      } else {
        print(paste("ERROR: Duplicate detected in sample", i, sep=" "))
      }
      
      rowNum <- rowNum + 1
    }
  }
}
