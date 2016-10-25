# diVaMo ("prova", "/Users/aingrosso/Dropbox/Progetti/IASI-CNR/Sviluppo/", "Input_Matrice_Campioni.csv", "Input_Matrice_Pathways.csv", "Classe", c("CASO", "CONTROLLO"), 1.1)
#diVaMo ("prova", "/Users/aingrosso/Dropbox/Progetti/IASI-CNR/Sviluppo/Matrici/", "TraspostaLPSMCAO.csv", , "Classe", c("Case", "Control"), 1.1)
#library(microbenchmark)

diVaMo <- function (name, path, sample_Matrix, pathways_Matrix, classId, classValues, threshold)
{
	#**************************************************************
	#START DEVELOPING PREPROCESSING FUNCTION
	#This function preprocesses samples and pathways, it will remove columns where genes' strenght is below the threshold.
	#**************************************************************
	preProcessing <- function(mC, stringPathology) {
		#print(labels(mC)[[2]][1])
		#print(subset(mC, (mC[,idClasse] == listaValoriClasse[1])))
		relationMatrix <- integer()
		arrayMaxHGenes <- apply(subset(mC, ((mC[,idClasse] == listaValoriClasse[1]) & (mC[,"Pathology"] == stringPathology))), 2, max)
		arrayMinHGenes <- apply(subset(mC, ((mC[,idClasse] == listaValoriClasse[1]) & (mC[,"Pathology"] == stringPathology))), 2, min)
		arrayMaxDGenes <- apply(subset(mC, ((mC[,idClasse] == listaValoriClasse[2]) & (mC[,"Pathology"] == stringPathology))), 2, max)
		arrayMinDGenes <- apply(subset(mC, ((mC[,idClasse] == listaValoriClasse[2]) & (mC[,"Pathology"] == stringPathology))), 2, min)
		for (i in 1:length(mC[1,])) {
			currGeneLabel <- colnames(mC)[i]
			if (currGeneLabel != "Pathology" && currGeneLabel != idClasse) {
				indexCurrGene <- grep(currGeneLabel, colnames(mC))
				if (
					!((as.numeric(arrayMaxDGenes[indexCurrGene])/as.numeric(arrayMinHGenes[indexCurrGene])) <= threshold &&
					(as.numeric(arrayMaxHGenes[indexCurrGene])/as.numeric(arrayMinDGenes[indexCurrGene])) <= threshold)
				) {
					#do something
					relationMatrix <- cbind(relationMatrix, c(i, grep(currGeneLabel, colnames(mPathways))[1]))
					colnames(relationMatrix)[ncol(relationMatrix)] <- currGeneLabel
				}
			}
			#cat("\r", as.integer((i/length(mC[1,])*100)), "%")
		}
		return(relationMatrix)
	}
	#**************************************************************
	# END DEVELOPING PREPROCESSING FUNCTION
	#**************************************************************
	
	#**************************************************************
	# START DEVELOPING LOG FUNCTION
	#**************************************************************
	insLog <- function (message) {
		if (newRow) {
			cat(currLogRow, format(Sys.time(), "%H:%M:%S"), message, sep=" ")
			cat(currLogRow, format(Sys.time(), "%H:%M:%S"), message, sep=" ", file = logFile, append = TRUE)
			currLogRow <<- currLogRow + 1
			newRow <<- FALSE
		}
		else {
			cat(message)
			cat(message, file = logFile, append=TRUE)
		}
		
		if (grepl("\n", c(message))) newRow <<- TRUE
	}
	#**************************************************************
	# END DEVELOPING LOG FUNCTION
	#**************************************************************

	#**************************************************************
	# START DEVELOPING EXIT FUNCTION
	#**************************************************************
	exitExe <- function() {
		close(logFile)
		cat(logStdMessage)
		return(FALSE)
	}
	#**************************************************************
	# END DEVELOPING EXIT FUNCTION
	#**************************************************************


#******************************************************************
#PROJECT NAME: DiVaMo
#VERSION: 0.6
#CREATION DATE: 22/03/2016
#LAST UPDATE: 13/07/2016
#AUTHOR: ALESSANDRO INGROSSO a.ingrosso@me.com
#******************************************************************
#Get OS and set path separator
#if (Sys.info()[['sysname']] == "Windows")
#	pathSep <- "\\"
#else
	pathSep <- "/"
	
#if path equal to "." or "" or blank, use the actual working dir
if (missing(path) || (path == ".") || (path == ""))
	path = getwd()
#Verifying that the path exists.
else if (file.exists(path))
	setwd(path)
else {
	cat("E01, cannot open path, no such directory\n")
	return(FALSE)
}

#append / to the path if absent
if (substr(path, nchar(path), nchar(path)) != pathSep)
	path = paste(path, pathSep, sep="")

#Initilizing log variables
currLogRow = 1
newRow = TRUE
logStdMessage = paste("A log file was created in ", path, name, "_log.txt\n", sep="")

#Assignments and checks
pathCampioni = sample_Matrix
if (missing(pathways_Matrix))
	pathPathways = NULL
else
	pathPathways = pathways_Matrix
idClasse = classId
listaValoriClasse <- classValues

#create log file
logFile <- file(paste(path, name, "_log.txt", sep=""), open = "wt")

insLog("Checking threshold parameter...")
if (threshold >= 1) {
	paramSoglia = threshold
} else {
	insLog("Error\n\tError 07: the threshhold parameter must be equal or greater than 1.\n")
	return(exitExe())
}
insLog("OK\n")

# Print all the single parameters value
insLog("\t\tINPUT PARAMETERS\n")
insLog(paste("\t\tElaboration name    : ", name,"\n"))
insLog(paste("\t\tWorking Directory   : ", path,"\n"))
insLog(paste("\t\tSamples matrix file : ", pathCampioni,"\n"))
if (is.null(pathPathways))
	insLog("\t\tPathways matrix file: empty\n")
else
	insLog(paste("\t\tPathways matrix file: ", pathPathways,"\n"))
insLog(paste("\t\tClass label name    : ", idClasse, "\n"))
insLog("\t\tClass values        : ")
for (i in 1:length(listaValoriClasse))
	insLog(paste(listaValoriClasse[i], " "))
insLog("\n")
insLog(paste("\t\tThreshold value     : ", paramSoglia, "\n"))

#******************************************************************
#INIZIO SVILUPPO REQUISITO R001
# mCampioni: matrice di campioni
# tableControlloR001: data frame che contiene l'unione dei valori delle colonne Pathology e Classe
#******************************************************************
#Caricamento matrice dei campioni
#print(table(count.fields(pathCampioni, sep=";")))

insLog("Reading samples matrix...")
#Verifying that samples matrix file exists
if (file.exists(pathCampioni)) {
	#Nota: è importante gestire diversi separatori e formati csv
	#Is the sample matrix malformed?
	mCampioni <- tryCatch(read.csv2(pathCampioni, row.names=1, header = TRUE),
						error = function(c) {
							insLog("Error\n\tError E03:")
							insLog(c$message)
						}
					)
	#Is the samples matrix empty?
	if (is.null(mCampioni)) return(exitExe())
} else {
	insLog("Error\n\tError E01: cannot open path, no such file.\n")
	return(exitExe())
}
insLog("OK\n")

#Matrix Summary
insLog("\t\tSAMPLES MATRIX SUMMARY\n")
insLog(paste("\t\tTotal samples:",length(mCampioni[,1]),"\n"))
for (i in 1:length(unique(mCampioni$Pathology))) {
	pathologyClasses <- unique(subset(mCampioni, Pathology == unique(mCampioni$Pathology)[i])[,idClasse])
	for (j in 1:length(pathologyClasses))
		insLog(paste("\t\tTotal samples labeled as",pathologyClasses[j],"for pathology", unique(mCampioni$Pathology)[i],":",length(subset(mCampioni, (Pathology == unique(mCampioni$Pathology)[i] & Classe == pathologyClasses[j]))[,idClasse]),"\n"))
}
insLog(paste("\t\tTotal gene attributes:",length(mCampioni[1,])-2,"\n")) #-2 because two matrix attributes are pathology and class

#R001: Unisco le colonne Pathology e Classe in un nuovo data frame, elimino i duplicati e verifico che il numero sia pari, cioè ad ogni patologia corrisponda almeno un caso e un controllo
insLog("Checking if every Pathologies have Case and Control samples...")
tableControlloR001 <- data.frame(paste(mCampioni$Pathology, mCampioni[,idClasse], sep="_"), row.names = rownames(mCampioni))

if (length(unique(duplicated(tableControlloR001[,1])))%%2==0)
	insLog("OK\n") 
else {
	insLog("Error\n\tError 04: missing at least one sample of the classe's values")
	return(exitExe())
}
#******************************************************************
#FINE SVILUPPO REQUISITO R001
#******************************************************************

#******************************************************************
#INIZIO SVILUPPO REQUISITO R002
# mPathways: matrice dei Pathways
# numPathology: lista delle patologie presenti nella matrice dei campioni
# pathways_Geni: matrice di output
# pathways_Score
#******************************************************************
#PASSO 1
#Is Pathways parameter filled?
if (is.null(pathPathways)) {
	insLog("Warning W02: pathways file not passed as parameter, all genes will be used.\n")
	#Extractin all genes labels from samples matrix
	subsetMSamples <- subset(labels(mCampioni)[[2]], (labels(mCampioni)[[2]] != "Pathology" & labels(mCampioni)[[2]] != idClasse))
	#Creating general pathway matrix
	mPathways = matrix(1,1,length(subsetMSamples), dimnames = list( c("General Pathway"), subsetMSamples ))
}
else {
	#Is Pathway matrix empty?
	mPathways = read.csv2(pathPathways, row.names=1, header = TRUE)
}

#SUMMARY PATHWAYS MATRIX
insLog("\t\tPATHWAYS MATRIX SUMMARY\n")
insLog(paste("\t\tTotal pathways:",length(mPathways[,1]),"\n"))
insLog(paste("\t\tTotal genes:",length(mPathways[1,]),"\n"))

#PASSO 2 R002
numPathology <- unique(mCampioni$Pathology)
#print(apply(expand.grid(labels(mPathways)[[1]],numPathology), 1, paste, collapse="_"))
pathways_Geni <- matrix(NA, nrow=(nrow(mPathways)*length(numPathology)), ncol=ncol(mPathways)+1, byrow=TRUE, dimnames=list(apply(expand.grid(labels(mPathways)[[1]],numPathology), 1, paste, collapse="_"), append(labels(mPathways)[[2]], "Pathology")))

#PASSO 3
pathways_Score <- matrix(NA, nrow=nrow(mPathways), ncol=length(numPathology), byrow=TRUE, dimnames=list(labels(mPathways)[[1]], numPathology))
pathways_Score <- cbind(pathways_Score, Size = NA)

#******************************************************************
#FINE SVILUPPO REQUISITO R002
#******************************************************************

#Preprocessing, reducing samples and pathways matrixes to the only valuable genes
sample_tbz_pathway <- list()
for (j in 1:length(numPathology)) {
	insLog(paste("Pre-processing", numPathology[j], "pathology in progress..."))
	sample_tbz_pathway[[j]] <- preProcessing(mCampioni, numPathology[j])
	insLog("OK\n")
	insLog(paste("The", as.integer(100-(ncol(sample_tbz_pathway[[j]])/length(mPathways[1,])*100)), "% of useless genes were removed by pre-processing activity.\n", ncol(sample_tbz_pathway[[j]]), "no discarded genes.\n"))
}

#******************************************************************
#INIZIO SVILUPPO REQUISITO R003
# pathwayGenes: vettore con etichette geni del pathway corrente
# keywordCaso: ad esempio "MCAO_CASO"
# keywordControllo: ad esempio "MCAO_CONTROLLO"
# nCasi: numero di CASI per patologia
# nCcontrolli: numero di CONTROLLI per patologia
# nAB: numero dei confronti tra campioni CASO e CONTROLLO
# contaUno: array con i conteggi degli 1 per ogni gene confrontato
# contaMenoUno: array con i conteggi dei -1 per ogni gene confrontato
# valEquilibrio: array con sbilanciamento
# valGood:
# matrice_Confronti_Good: contiene i confronti tra i valori di matrice_Confronti e i valGood per ogni gene
# valContaGoodUno: è un array che conta gli 1 presenti nelle singole righe della matrice_Confronti_Good
# valContaGoodMenoUno: è un array che conta i -1 presenti nelle singole righe della matrice_Confronti_Good
# minGoodUno: è il minimo valore tra tutti i valori dell'array valContaGoodUno
# minGoodMenoUno: è il minimo valore tra tutti i valori dell'array valContaGoodMenoUno
# contaGoodZero: è il numero di zeri presenti nell'array valContaGoodUno
# contaGoodMenoZero: è il numero di zeri presenti nell'array valContaGoodMenoUno
# maxEquilibrio: è il massimo valore tra tutti i valori contenuti nell'array valEquilibrio
# alpha: è il massimo tra minGoodUno e minGoodMenoUno
#******************************************************************
for(i in 1:nrow(mPathways)) { #PASSO 1
	
	#PASSO 2
	for (j in 1:length(numPathology)) {
		pathwayGenes <- integer()
		#print(labels(mPathways)[[1]][i]) #stampa label riga
		for (k in 1:ncol(mPathways)) { #PASSO 2
			if (mPathways[i,k]==1){
				#print(colnames(mPathways)[k])
				#print(colnames(sample_tbz_pathway[[1]]))
				#print(grep(colnames(mPathways)[k], colnames(sample_tbz_pathway[[j]]), fixed = TRUE)[1])
				#Extract index position of the current gene from sample_tbz_pathway 
				indexColSampleTbzPathway = grep(colnames(mPathways)[k], colnames(sample_tbz_pathway[[j]]), fixed = TRUE)[1]
				#print(sample_tbz_pathway[[1]][1,indexColSampleTbzPathway]) #estrae posizione nella matrice campioni
				#if index is not NA, it means that the gene was non discarded during preprocessing, add index in the pathwayGenes
				if (!is.na(indexColSampleTbzPathway)) {
					pathwayGenes <- append(pathwayGenes, indexColSampleTbzPathway)
				}
			}
		}
		
		if (length(pathwayGenes)!=0) {
		
			insLog(paste("Building", row.names(mPathways)[i], "pathway,", numPathology[j], "pathology\n", sep=" "))
			
			#PASSO 3
			keywordCaso = paste(numPathology[j], listaValoriClasse[1], sep="_")
			keywordControllo = paste(numPathology[j], listaValoriClasse[2], sep="_")
			nCasi = table(tableControlloR001 == keywordCaso)["TRUE"]
			nControlli = table(tableControlloR001 == keywordControllo)["TRUE"]
			nAB = nCasi * nControlli
			insLog(paste("Numero totale confronti",nAB,"\n"))
			matrice_Confronti <- array(NA, c(nAB,length(pathwayGenes)))
			rigaCorrenteConfronto = 1
			
			#PASSO 4
			#complessità n^3
			#for (k in match(keywordCaso,tableControlloR001):(nCasi-1+match(keywordCaso,tableControlloR001))) {
			for (k in 1:nrow(tableControlloR001)) {
				if (keywordCaso == tableControlloR001[k,1]) {
					insLog(paste("Comparing Case", rownames(tableControlloR001)[[k]], "with Controls:\n"))
					#for (n in match(keywordControllo,tableControlloR001):(nControlli-1+match(keywordControllo,tableControlloR001))) {
					for (n in 1:nrow(tableControlloR001)) {
						if (keywordControllo == tableControlloR001[n,1]) {
							insLog(paste("\t", rownames(tableControlloR001)[[n]], "\n"))
							#options(digits.secs=6)
							for (m in 1:length(pathwayGenes)) {
								#insLog(paste("\t\tGene:", pathwayGenes[m]))
								#t1 = Sys.time()
								if (!is.null(mCampioni[k,sample_tbz_pathway[[j]][1,pathwayGenes[m]]])) {
									if ((mCampioni[k,sample_tbz_pathway[[j]][1,pathwayGenes[m]]]/mCampioni[n,sample_tbz_pathway[[j]][1,pathwayGenes[m]]]) > paramSoglia)
										matrice_Confronti[rigaCorrenteConfronto,m] = 1
									else if ((mCampioni[n,sample_tbz_pathway[[j]][1,pathwayGenes[m]]]/mCampioni[k,sample_tbz_pathway[[j]][1,pathwayGenes[m]]]) > paramSoglia)
										matrice_Confronti[rigaCorrenteConfronto,m] = -1
									else
										matrice_Confronti[rigaCorrenteConfronto,m] = 0
								}
								#t2 = Sys.time()
								#cat("Durata",(t2-t1),"\n")
							}
							rigaCorrenteConfronto = rigaCorrenteConfronto + 1
						}
					}
				}
			}
			#return(exitExe())
			
			contaUno <- integer()
			contaMenoUno <- integer()
			valEquilibrio <- integer()
			valGood <- integer()
			valContaGoodUno <- integer()
			valContaGoodMenoUno <- integer()
			for (k in 1:length(pathwayGenes)) {
				#PASSO 5
				if (is.na(table(matrice_Confronti[,k] == 1)["TRUE"]))
					contaUno <- append(contaUno, 0)
				else
					contaUno <- append(contaUno, table(matrice_Confronti[,k] == 1)["TRUE"])
				
				if (is.na(table(matrice_Confronti[,k] == -1)["TRUE"]))
					contaMenoUno <- append(contaMenoUno, 0)
				else {
					contaMenoUno <- append(contaMenoUno, table(matrice_Confronti[,k] == -1)["TRUE"])
					#sum = sum + table(matrice_Confronti[,k] == -1)["TRUE"]
				}
				
				#PASSO 6
				if (contaUno[k] == 0)
					valEquilibrio <- append(valEquilibrio, contaMenoUno[k])
				else if (contaMenoUno[k] == 0)
					valEquilibrio <- append(valEquilibrio, contaUno[k])
				else
					valEquilibrio <- append(valEquilibrio, 0)
				
				#PASSO 7
				if ((contaUno[k] == 0) & (contaMenoUno[k] > 0))
					valGood = append(valGood, -1)
				else if ((contaUno[k] > 0) & (contaMenoUno[k] == 0))
					valGood = append(valGood, 1)
				else
					valGood = append(valGood, 0)
			}
			#print(paste(contaUno, pathwayGenes, sep="_"))
			#print(paste(contaMenoUno, pathwayGenes, sep="_"))
			#print(sum)
			#print(valEquilibrio)
			#print(valGood)
			
			#PASSO 8
			matrice_Confronti_Good <- array(NA, c(nAB,length(pathwayGenes)))
			for (k in 1:nrow(matrice_Confronti)) {
				for (n in 1:ncol(matrice_Confronti)) {
					if (!is.na(matrice_Confronti[k,n])) {
						if ((matrice_Confronti[k,n] * valGood[n]) > 0)
							matrice_Confronti_Good[k,n] = 1
						else
							matrice_Confronti_Good[k,n] = 0
					}
				}
				if (!is.na(table(matrice_Confronti_Good[k,] == 1)["TRUE"]))
					valContaGoodUno = append(valContaGoodUno, table(matrice_Confronti_Good[k,] == 1)["TRUE"])
				else
					valContaGoodUno = append(valContaGoodUno, 0)
				if (!is.na(table(matrice_Confronti_Good[k,] == -1)["TRUE"]))
					valContaGoodMenoUno = append(valContaGoodMenoUno, table(matrice_Confronti_Good[k,] == -1)["TRUE"])
				else
					valContaGoodMenoUno = append(valContaGoodMenoUno, 0)
			}
			#print(matrice_Confronti_Good)
			#print(valContaGoodUno)
			#print(valContaGoodMenoUno)
			
			#PASSO 9
			minGoodUno = min(valContaGoodUno)
			minGoodMenoUno = min(valContaGoodMenoUno)
			contaGoodZero = table(valContaGoodUno == 0)["TRUE"]
			if (is.na(contaGoodZero)) contaGoodZero = 0
			contaGoodMenoZero = table(valContaGoodMenoUno == 0)["TRUE"]
			if (is.na(contaGoodMenoZero)) contaGoodMenoZero = 0
			maxEquilibrio = max(valEquilibrio)
			alpha = max(minGoodUno,minGoodMenoUno)
			#print(minGoodUno)
			#print(minGoodMenoUno)
			#print(contaGoodZero)
			#print(contaGoodMenoZero)
			#print(maxEquilibrio)
			#print(alpha)
			
			#PASSO 10
			if (minGoodUno > minGoodMenoUno)
				nalpha = contaGoodZero
			else if (minGoodMenoUno > minGoodUno)
				nalpha = contaGoodMenoZero
			else
				nalpha = min(contaGoodZero,contaGoodMenoZero)
			
			#PASSO 11
			for (k in 1:length(pathwayGenes)) {
				pathways_Geni[paste(labels(mPathways)[[1]][i],numPathology[j],sep="_"), sample_tbz_pathway[[j]][2,pathwayGenes[k]]] = (max(abs(contaUno[k]),contaMenoUno[k]))/nAB
			}
			pathways_Geni[paste(labels(mPathways)[[1]][i],numPathology[j],sep="_"),"Pathology"] = numPathology[j]
			#print(paste(labels(mPathways)[[1]][i],numPathology[j],sep="_"))
			#print(pathways_Geni[paste(labels(mPathways)[[1]][i],numPathology[j],sep="_"),])
			
			#PASSO 12
			pathways_Score[i,j] = ((alpha * nAB) + (nAB - nalpha))/nAB
			pathways_Score[i,"Size"] = length(pathwayGenes)
		
		}
	
	}
	#print(pathways_Geni)
	#print(paste(nalpha, contaGoodZero, "\n"))
}
#SEGUE PASSO 11
#elaboro totali per patologia
#Doing sums per columns only if the matrix has more than one row
for (i in 1:length(numPathology)) {
	#Doing sum per columns only if there is more than one row per pathology
	if (length(pathways_Geni[pathways_Geni[,"Pathology"] == i,"Pathology"]) > 1) {
		pathways_Geni <- rbind(pathways_Geni, Totale = colSums(pathways_Geni[pathways_Geni[,"Pathology"] == i,], na.rm=TRUE))
		pathways_Geni[nrow(pathways_Geni),"Pathology"] = i
	}
}
	
#ordino le righe per patologia
pathways_Geni[order(pathways_Geni[,"Pathology"]),]
write.csv2(t(pathways_Geni), paste(path, name, "_Pathways_Genes_Matrix.csv", sep=""))#"/Users/aingrosso/Dropbox/Progetti/IASI-CNR/Sviluppo/Matrice_Pathways_Geni.csv")
write.csv2(pathways_Score, paste(path, name, "_Pathways_Score_Matrix.csv", sep=""))#"/Users/aingrosso/Dropbox/Progetti/IASI-CNR/Sviluppo/Matrice_Pathways_Score.csv")

#******************************************************************
#FINE SVILUPPO REQUISITO R003
#******************************************************************
insLog(paste("Elaboration completed\n", logStdMessage))
close(logFile)
return(TRUE)
}