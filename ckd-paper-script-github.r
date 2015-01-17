# clinical-cohort-compare.r
# code author: Craig Smail
#
# description: 
# 1./ construct a standardized matrix of clinical lab values and perform hierarchical clustering
# of resulting matrix to find cohorts of interest;
# 2./ construct a matched control/comparison cohort using cohort from (1.) as reference - matched on
# BMI, gender, and age;
# 3./ compute mean values, significance tests, and rates of change for each clinical lab for each group

# set working directory
# parameters: directory containing input data
setwd("<<dir>>")

# import clinical data
# parameters: input files (tab separated, with headers, linked with unique person ID across all files)
# minimum column requirements (case sensitive):
# - findingTable: IDPerson, FindValNum, IDValue, FindDate
# - personTable: IDPerson, DOB (format: YYYY-MM-DD 00:00:00.000), Gender (M/F)
# - labsTable: IDPerson, IDValue, ResultValueNum
# - dxTable: IDPerson, DxOnsetDate
findingTable = read.csv("findings.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE) # blood pressure, height, weight
personTable = read.csv("person.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE) # year of birth, gender
labsTable = read.csv("labs.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE) # eGFR, creatinine, etc.
dxTable = read.csv("dx.txt", sep="\t", header=TRUE, stringsAsFactors=FALSE) # ICD-9 code (or equivalent), date of dx onset

# search for cases with disease of interest
cases = grep("<<ICD9 code>>",dxTable$DxCode) # find patients with ICD-9 of interest
casesID = as.vector(dxTable$IDPerson[cases]) # extract patient ID
casesProfile = dxTable[which(dxTable$IDPerson %in% casesID), ] # extract all dx for case patients
casesLimit = casesProfile[grep("<<ICD9 code(s)", casesProfile$DxCode), ] # filter only to dx(s) of interest
casesFilter = casesLimit[which(casesLimit$IDPerson %in% casesID), ] 

# get all unique cases and assign unique ID
uniquePatients = unique(casesFilter$IDPerson) 
uniquePatientsLength = length(unique(casesFilter$IDPerson)) 
casesFilter$UniqueID = NA
for (i in 1:uniquePatientsLength) {
    tempMatch = which(casesFilter$IDPerson==uniquePatients[i])
	casesFilter$UniqueID[tempMatch] = i 
}

# create labs matrix (all observations in labs)
# parameters: lab names of interest from input data
# optional: use all unique labs in dataset
uniquePat = unique(labsTable$IDPerson)
uniqueLabName = c(202, 2) # lab names of interest
# uniqueLabName = unique(labsTable$IDValue)

# get counts of each unique lab name
collector = array(NA,c(length(uniqueLabName), 2, length(uniquePat)))
for (i in 1:length(uniquePat)) {
  for (j in 1:length(uniqueLabName)) {
    collector[j, 1, i] = uniqueLabName[j]
    collector[j, 2, i] = length(which(labsTable$IDValue==uniqueLabName[j] & labsTable$IDPerson==uniquePat[i]))
  }
}

# determine sparseness of data for each lab test in dataset
collectorSparse = array(NA,c(length(uniqueLabName), 3))
for (i in 1:length(uniqueLabName)) {
	collectorSparse[i, 1] = uniqueLabName[i]
	collectorSparse[i, 2] = sum(as.numeric(collector[i, 2, ]))
	collectorSparse[i, 3] = length(which(as.numeric(collector[i, 2, ])==0))
}

# Data standardization
# parameters: labs passing user-defined completeness (e.g. use labs with >= 1 count across all patients)
thresholdLabName = c(<<lab string(s), comma separated>>) 
collectorLabPop = array(NA,c(length(uniquePat), length(thresholdLabName) + 1))

# get count for each patient, for each lab test that meets threshold 
for (i in 1:length(uniquePat)) {
	for (j in 1:length(thresholdLabName)) {
		matchLab = length(which(labsTable$IDPerson==uniquePat[i] & labsTable$IDValue==thresholdLabName[j]))
		collectorLabPop[i, 1] = uniquePat[i]
		offsetIndex = j+1
		collectorLabPop[i, offsetIndex] = ifelse(matchLab > 0, matchLab, NA) 
	}
}

uniquePatients1 = unique(collectorLabPop[complete.cases(collectorLabPop), 1]) # vector of unique patient IDs with > 0 lab results for each lab test

# get mean lab values for each lab test, for each patient, and add mean blood test (diastolic and systolic) and age in years (2014 - year of birth)
# parameters: 
# - BPS (blood pressure systolic IDValue)
# - BPD (blood pressure diastolic IDValue)
labsMatrix = matrix(nrow=length(uniquePatients1), ncol=length(thresholdLabName) + 5, byrow=TRUE)
colnames(labsMatrix) = c("ID", "StandID", "Age", thresholdLabName, "BPS", "BPD") # dim column names
labsMatrix[ ,1] = uniquePatients1
labsMatrix[ ,2] = seq(from=1, to=length(uniquePatients1), by=1) # standardized patient ID, to reduce from eight digit ID found in dataset

today = Sys.Date()
for (i in 1:length(uniquePatients1)) {
	for (j in 1:length(thresholdLabName)) {
		labValueMean = mean(as.numeric(labsTable$ResultValueNum[which(labsTable$IDPerson==uniquePatients1[i] & labsTable$IDValue==thresholdLabName[j])]))
		offsetIndex2 = j+3
		labsMatrix[i, offsetIndex2] = labValueMean
		labsMatrix[i, 3] = as.numeric(today - as.Date(personTable$PsnDOB[match(uniquePatients1[i], personTable$IDPerson)])) / 365  # get YOB and compute approx age in years
		bps = mean(as.numeric(findingTable$FindValNum[which(findingTable$IDPerson==uniquePatients1[i] & findingTable$IDValue==<<BPS string>>)])) # BPS 
		bpd = mean(as.numeric(findingTable$FindValNum[which(findingTable$IDPerson==uniquePatients1[i] & findingTable$IDValue==<<BPD string>>)])) # BPD
		labsMatrix[i, ncol(labsMatrix) - 1] = ifelse(is.null(bps)==TRUE, NA, bps) # if no BPS value, write NA
		labsMatrix[i, ncol(labsMatrix)] = ifelse(is.null(bpd)==TRUE, NA, bpd) # if no BPD value, write NA
	}
}

labsMatrix = labsMatrix[complete.cases(labsMatrix), ] # remove patients with NaN

# standardize variables by calculating z-score
# calculate absolute deviation
sfArray = array(NA, c(1, ncol(labsMatrix)-3))

meanDevCollect = list()
for (i in 1:length(sfArray)) {
	meanDevCollect[[i]] = 0
}

for (j in 4:ncol(labsMatrix)) {
	k = j - 3
	meanDev = labsMatrix[, j] - mean(labsMatrix[, j])
	meanDevCollect[[k]] = abs(meanDev)
	sfArray[k] = (1/nrow(labsMatrix)) * sum(meanDevCollect[[k]])
}

# calculate Z-score
labsMatrixStand = labsMatrix[ , -c(2,3)] # remove age and standardized ID columns (not appropriate for inclusion in clustering)
for (j in 1:nrow(labsMatrixStand)) {
	for (i in 2:ncol(labsMatrixStand)) {
		labsMatrixStand[j, i] = (labsMatrixStand[j, i] - mean(labsMatrix[, i+2])) / sfArray[i-1]
	}
}

# compute distance and fit clustering model 
# parameters: change distance and linkage methods to suit application
d = dist(labsMatrixStand[, -1], method="euclidean") # remove patient ID for cluster
hCluster = hclust(d, method="complete")
dend = as.dendrogram(hCluster)

# heat map
rowCol = array(NA,c(nrow(labsMatrixStand)))

for (i in 1:nrow(labsMatrixStand)) {
	rowCol[i] = ifelse(labsMatrixStand[i, 1] %in% casesID, "#ff69b4", "#0000FF") # highlight cases
}

labNamesClean = c(thresholdLabName, "BPS", "BPD") # for figure legend

# red = lower values, yellow = higher values
pdf("heatMap.pdf", width=40, height=15)
heatmap(labsMatrixStand[, -1], Rowv=dend, Colv=NULL, ColSideColors=rainbow(length(labNamesClean), start=0.3, end=1), keep.dendro=FALSE,
			RowSideColors=rowCol, labCol=labNamesClean, labRow=labsMatrix[, 2], cexRow=0.1)
dev.off()

clusterFocus = as.vector(read.delim("<<observations of interest>>", header=FALSE)) # save vector of patient IDs from cluster(s) of interest (from manual inspection)
patientFocus = labsMatrix[which(labsMatrix[ ,2] %in% unlist(clusterFocus)), ] # filter lab results matrix to include only those patients in the cluster of interest
colnames(patientFocus) = c("ID", "StandID", "Age", thresholdLabName, "BPS", "BPD") # dim column names
patientCompare = labsMatrix[which(!labsMatrix[ ,2] %in% unlist(clusterFocus)), ] # filter lab results matrix to include every patient not in cluster of interest

# Summary lab statistics
## Case Group demographics, labs, and dx values

patientsN = nrow(patientFocus)

# summaryStats function
# prints gender count in cohort, and lab means and standard deviation
# parameters: 
# cohort (input matrix containing lab data for cohort of interest)
summaryStats = function(cohort) {
	print(c("gender male: ", nrow(personTable[which(personTable$IDPerson %in% cohort[ ,1] & personTable$PsnGender=="M"),])))
	print(c("gender female: ", nrow(personTable[which(personTable$IDPerson %in% cohort[ ,1] & personTable$PsnGender=="F"),])))
	for (i in 3:ncol(cohort)) {
		print(c(colnames(cohort)[i], " mean: ", mean(labsMatrix[which(labsMatrix[ ,1] %in% cohort[, 1]), i])))
		print(c(colnames(cohort)[i], " sd: ", sd(labsMatrix[which(labsMatrix[ ,1] %in% cohort[, 1]), i])))
	}
}

summaryStats(patientFocus)

# calculate BMI (imperial formula: (weight in pounds * 703)/(height in inches)^2)
# parameters:
# - height (string representing height measurements in dataset)
# - weight (string representing weight measurements in dataset)
bmiCollector = array(NA, c(nrow(labsMatrix), 2))
bmiCollector[, 1] = labsMatrix[ ,1]

for (i in 1:nrow(bmiCollector)) {
	height = findingTable[which(findingTable$IDPerson==bmiCollector[i, 1] & findingTable$IDValue==<<height string>>), "FindValNum"]
	weightAll = findingTable[which(findingTable$IDPerson==bmiCollector[i, 1] & findingTable$IDValue==<<weight string>>), c("FindValNum", "FindDate")]
	weightRecent = weightAll$FindValNum[which.max(as.Date(weightAll[, 2]))]
	bmiCollector[i, 2] = ifelse(length(height) > 0 & length(weightRecent) > 0, (weightRecent * 703) / mean(height)^2, NA) # get mean height to account for repeated measures
}

# Control group
# Match

caseControlMatch = list() # dim list for matching controls for each case

for (i in 1:length(patientFocus[ ,1])) {
 caseControlMatch[[i]] = 0
}

for (i in 1:nrow(patientFocus)) {
  caseGender = personTable$PsnGender[personTable$IDPerson==patientFocus[i, 1]] 
  caseBMI = bmiCollector[which(bmiCollector[, 1]==patientFocus[i, 1]), 2]
  controlMatchGender = personTable[which(personTable$IDPerson %in% patientCompare[ ,1] & personTable$PsnGender==caseGender), 1] 
  controlMatchBMI = bmiCollector[which(bmiCollector[ ,1] %in% patientCompare[ ,1] 
										& bmiCollector[, 2] < caseBMI + 8 &  bmiCollector[, 2] > caseBMI - 8)] 
  caseControlMatch[[i]] = controlMatchGender[which(controlMatchGender %in% controlMatchBMI)]
}

caseControlMatch2 = list() # dim list for matching controls for each case

for (i in 1:length(patientVector)) {
	caseControlMatch2[[i]] = sample(caseControlMatch[[i]], size=min(length(caseControlMatch[[i]]), 2)) # randomly sample no more than 2 matching controls for each case
}

controlMatchVector = as.vector(unlist(caseControlMatch2))
controlMatchVector = controlMatchVector[complete.cases(controlMatchVector)] 
controlMatchVector = unique(controlMatchVector) # check for duplicated patient IDs (i.e. matched control patients similar to > 1 case patient)

patientsNControl = length(controlMatchVector) 

# control summary stats
summaryStats(labsMatrix[which(labsMatrix[ ,1] %in% controlMatchVector), ])

# significance tests (t-tests, two-sided)
# prints p-value only
for (i in 3:ncol(labsMatrix)) {
	print(c(colnames(labsMatrix)[i], t.test(patientFocus[ ,3], labsMatrix[which(labsMatrix[ ,1] %in% controlMatchVector) , i], alternative="two.sided")$p.value))
}

# estimated annual rate of change

collectorRateChange <- array(NA,c(length(thresholdLabName), 5))
colnames(collectorRateChange) <- c("LabTest", "MeanRateChangeCase", "MeanRateChangeControl", "nCase", "nControl")

# case

rateChange = list()
for (i in 1:length(patientVector)) {
	rateChange[[i]] = 0
}

for (j in 1:length(thresholdLabName)) {
	for (i in 1:length(patientVector)) {
		currentPat = labsTable[which(labsTable$IDPerson==patientVector[i] & labsTable$IDValue==thresholdLabName[j]), c("ResultValueNum", "ResultDate")]
		oldLab = as.numeric(currentPat[which.min(as.Date(currentPat[, 2])), 1])
		newLab = as.numeric(currentPat[which.max(as.Date(currentPat[, 2])), 1])
		dateDif = as.Date(currentPat[which.max(as.Date(currentPat[, 2])), 2])  - 
			as.Date(currentPat[which.min(as.Date(currentPat[, 2])), 2])
		rate = (((newLab - oldLab) / unlist(as.numeric(dateDif))) * 365) # calculated annual rate of change
		rateChange[[i]] = ifelse(is.nan(rate), 0, rate)
	}

	meanRateChange = unlist(rateChange)
	collectorRateChange[j, 1] = labClean[j]
	collectorRateChange[j, 2] = mean(meanRateChange)
	collectorRateChange[j, 4] = length(meanRateChange)
}

# control

rateChangeControl = list()
for (i in 1:length(controlMatchVector)) {
	rateChangeControl[[i]] = 0
}

for (j in 1:length(thresholdLabName)) {
	for (i in 1:length(controlMatchVector)) {
		currentPat = labsTable[which(labsTable$IDPerson==controlMatchVector[i] & labsTable$IDValue==thresholdLabName[j]), c("ResultValueNum", "ResultDate")]
		oldLab = as.numeric(currentPat[which.min(as.Date(currentPat[, 2])), 1])
		newLab = as.numeric(currentPat[which.max(as.Date(currentPat[, 2])), 1])
		dateDif = as.Date(currentPat[which.max(as.Date(currentPat[, 2])), 2])  - 
			as.Date(currentPat[which.min(as.Date(currentPat[, 2])), 2])
		rate = (((newLab - oldLab) / unlist(as.numeric(dateDif))) * 365) # calculated annual rate of change
		rateChangeControl[[i]] = ifelse(is.nan(rate), 0, rate)
	}

	meanRateChange = unlist(rateChangeControl)
	collectorRateChange[j, 1] = labClean[j]
	collectorRateChange[j, 3] = mean(meanRateChange)
	collectorRateChange[j, 5] = length(meanRateChange)
}