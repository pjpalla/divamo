\name{diVaMo}
\alias{diVaMo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
DiVaMo
}
\description{
  This library compares two specimens' classes, desease and healty, and it uses pathways to collect a specific genes' subset.
  Desease specimens are named case, healty ones are named control.
}
\usage{
diVaMo()
}
%- maybe also 'usage' for other objects documented here.
\details{
  No more details are needed.
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{comp1 }{Description of 'comp1'}
  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
 ~put references to the literature/web site here ~
}
\author{
  Alessandro Ingrosso
}
\note{
  No further notes needed.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
 ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function () 
{
    pathCampioni = "/Users/aingrosso/Dropbox/Progetti/IASI-CNR/Sviluppo/Input_Matrice_Campioni.csv"
    pathPathways = "/Users/aingrosso/Dropbox/Progetti/IASI-CNR/Sviluppo/Input_Matrice_Pathways.csv"
    idClasse = "Classe"
    listaValoriClasse <- c("CASO", "CONTROLLO")
    paramSoglia = 1.1
    mCampioni = read.csv2(pathCampioni, row.names = 1, header = TRUE)
    tableControlloR001 <- paste(mCampioni$Pathology, mCampioni$Classe, 
        sep = "_")
    if (length(unique(tableControlloR001[duplicated(tableControlloR001)]))\%\%2 == 
        0) 
        print("controllo superato")
    else print("controllo non superato")
    mPathways = read.csv2(pathPathways, row.names = 1, header = TRUE)
    numPathology <- unique(mCampioni$Pathology)
    pathways_Geni <- matrix(NA, nrow = (nrow(mPathways) * length(numPathology)), 
        ncol = ncol(mPathways) + 1, byrow = TRUE, dimnames = list(paste(labels(mPathways)[[1]], 
            numPathology, sep = "_"), append(labels(mPathways)[[2]], 
            "Pathology")))
    pathways_Score <- matrix(NA, nrow = nrow(mPathways), ncol = length(numPathology), 
        byrow = TRUE, dimnames = list(labels(mPathways)[[1]], 
            numPathology))
    pathways_Score <- cbind(pathways_Score, Size = NA)
    for (i in 1:nrow(mPathways)) {
        pathwayGenes <- character()
        for (j in 1:ncol(mPathways)) {
            if (mPathways[i, j] == 1) 
                pathwayGenes <- append(pathwayGenes, c(labels(mPathways)[[2]][j]))
        }
        for (j in 1:length(numPathology)) {
            print(paste("Elaboro il Pathway", row.names(mPathways)[i], 
                "della patologia", numPathology[j], sep = " "))
            keywordCaso = paste(numPathology[j], listaValoriClasse[1], 
                sep = "_")
            keywordControllo = paste(numPathology[j], listaValoriClasse[2], 
                sep = "_")
            nCasi = table(tableControlloR001 == keywordCaso)["TRUE"]
            nControlli = table(tableControlloR001 == keywordControllo)["TRUE"]
            nAB = nCasi * nControlli
            matrice_Confronti <- array(NA, c(nAB, length(pathwayGenes)))
            rigaCorrenteConfronto = 1
            for (k in match(keywordCaso, tableControlloR001):(nCasi - 
                1 + match(keywordCaso, tableControlloR001))) {
                for (n in match(keywordControllo, tableControlloR001):(nControlli - 
                  1 + match(keywordControllo, tableControlloR001))) {
                  for (m in 1:length(pathwayGenes)) {
                    if (!is.null(mCampioni[k, pathwayGenes[m]])) {
                      if ((mCampioni[k, pathwayGenes[m]]/mCampioni[n, 
                        pathwayGenes[m]]) > paramSoglia) 
                        matrice_Confronti[rigaCorrenteConfronto, 
                          m] = 1
                      else if ((mCampioni[n, pathwayGenes[m]]/mCampioni[k, 
                        pathwayGenes[m]]) > paramSoglia) 
                        matrice_Confronti[rigaCorrenteConfronto, 
                          m] = -1
                      else matrice_Confronti[rigaCorrenteConfronto, 
                        m] = 0
                    }
                  }
                  rigaCorrenteConfronto = rigaCorrenteConfronto + 
                    1
                }
            }
            contaUno <- integer()
            contaMenoUno <- integer()
            valEquilibrio <- integer()
            valGood <- integer()
            valContaGoodUno <- integer()
            valContaGoodMenoUno <- integer()
            for (k in 1:length(pathwayGenes)) {
                if (is.na(table(matrice_Confronti[, k] == 1)["TRUE"])) 
                  contaUno <- append(contaUno, 0)
                else contaUno <- append(contaUno, table(matrice_Confronti[, 
                  k] == 1)["TRUE"])
                if (is.na(table(matrice_Confronti[, k] == -1)["TRUE"])) 
                  contaMenoUno <- append(contaMenoUno, 0)
                else {
                  contaMenoUno <- append(contaMenoUno, table(matrice_Confronti[, 
                    k] == -1)["TRUE"])
                }
                if (contaUno[k] == 0) 
                  valEquilibrio <- append(valEquilibrio, contaMenoUno[k])
                else if (contaMenoUno[k] == 0) 
                  valEquilibrio <- append(valEquilibrio, contaUno[k])
                else valEquilibrio <- append(valEquilibrio, 0)
                if ((contaUno[k] == 0) & (contaMenoUno[k] > 0)) 
                  valGood = append(valGood, -1)
                else if ((contaUno[k] > 0) & (contaMenoUno[k] == 
                  0)) 
                  valGood = append(valGood, 1)
                else valGood = append(valGood, 0)
            }
            matrice_Confronti_Good <- array(NA, c(nAB, length(pathwayGenes)))
            for (k in 1:nrow(matrice_Confronti)) {
                for (n in 1:ncol(matrice_Confronti)) {
                  if (!is.na(matrice_Confronti[k, n])) {
                    if ((matrice_Confronti[k, n] * valGood[n]) > 
                      0) 
                      matrice_Confronti_Good[k, n] = 1
                    else matrice_Confronti_Good[k, n] = 0
                  }
                }
                if (!is.na(table(matrice_Confronti_Good[k, ] == 
                  1)["TRUE"])) 
                  valContaGoodUno = append(valContaGoodUno, table(matrice_Confronti_Good[k, 
                    ] == 1)["TRUE"])
                else valContaGoodUno = append(valContaGoodUno, 
                  0)
                if (!is.na(table(matrice_Confronti_Good[k, ] == 
                  -1)["TRUE"])) 
                  valContaGoodMenoUno = append(valContaGoodMenoUno, 
                    table(matrice_Confronti_Good[k, ] == -1)["TRUE"])
                else valContaGoodMenoUno = append(valContaGoodMenoUno, 
                  0)
            }
            minGoodUno = min(valContaGoodUno)
            minGoodMenoUno = min(valContaGoodMenoUno)
            contaGoodZero = table(valContaGoodUno == 0)["TRUE"]
            contaGoodMenoZero = table(valContaGoodMenoUno == 
                0)["TRUE"]
            maxEquilibrio = max(valEquilibrio)
            alpha = max(minGoodUno, minGoodMenoUno)
            if (minGoodUno > minGoodMenoUno) 
                nalpha = contaGoodZero
            else if (minGoodMenoUno > minGoodUno) 
                nalpha = contaGoodMenoZero
            else nalpha = min(contaGoodZero, contaGoodMenoZero)
            for (k in 1:length(pathwayGenes)) {
                pathways_Geni[paste(labels(mPathways)[[1]][i], 
                  numPathology[j], sep = "_"), pathwayGenes[k]] = (max(abs(contaUno[k]), 
                  contaMenoUno[k]))/nAB
            }
            pathways_Geni[paste(labels(mPathways)[[1]][i], numPathology[j], 
                sep = "_"), "Pathology"] = numPathology[j]
            pathways_Score[i, j] = ((alpha * nAB) + (nAB - nalpha))/nAB
            pathways_Score[i, "Size"] = length(pathwayGenes)
        }
    }
    for (i in 1:length(numPathology)) {
        pathways_Geni <- rbind(pathways_Geni, Totale = colSums(pathways_Geni[pathways_Geni[, 
            "Pathology"] == i, ], na.rm = TRUE))
        pathways_Geni[nrow(pathways_Geni), "Pathology"] = numPathology[j]
    }
    pathways_Geni[order(pathways_Geni[, "Pathology"]), ]
    write.csv2(pathways_Geni, "/Users/aingrosso/Dropbox/Progetti/IASI-CNR/Sviluppo/Matrice_Pathways_Geni.csv")
    write.csv2(pathways_Score, "/Users/aingrosso/Dropbox/Progetti/IASI-CNR/Sviluppo/Matrice_Pathways_Score.csv")
    print("finito")
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
