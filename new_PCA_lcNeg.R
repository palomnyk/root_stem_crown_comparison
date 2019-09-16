# New PCA plot

setwd(file.path('~','git','root_stem_crown_comparison'))

ds = 'lcNeg'

#read in meta-data
metadata = read.table(file.path('.','data','meta_data_no_whole.csv'), 
                      sep=",", 
                      header=TRUE, 
                      row.names = 1, 
                      check.names = FALSE,
                      stringsAsFactors=FALSE)

metabolites = read.table(file.path('.', 'data','lcNegNoWhole.csv'),
                         sep=",",
                         header=TRUE,
                         row.names = 1,
                         check.names = FALSE)

metabolites = log2(metabolites + 1)

#create PCA
group_pca = prcomp(t(na.omit(metabolites)), 
                   center = TRUE,
                   scale = TRUE)

#extract PCA matrix and convert to dataframe
myPCA = data.frame(group_pca$x)

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

#set color palette
palette( brewer.pal(7,"Accent") )

pdf('newPCA_lcNeg.pdf')
par(bty = 'l', 
    mar = c(5, 4, 5, 2) + 0.1,
    cex.lab = 1.25
)
plot(myPCA$PC1, myPCA$PC2,
     col='black',
     pch = c(24,21,22)[as.numeric(as.factor(metadata['structural location',]))],
     bg = as.factor(metadata['root vs stem',]),
     cex = 3,
     xlab = 'PCA 1, 31.71%',
     ylab = 'PCA 2, 22.32%',
     main = 'PCA of LC/MS Negative data'
)
# legend('topright', 
#        legend = c('crown root', 'lateral root', 'stem'),
#        pt.bg = as.factor(c('root', 'root', 'stem')),
#        pch = c(24,21,22),
#        cex = 1
# )
dev.off()


lmeMetadata = function(testingData, label){
  if(!require(nlme)){install.packages("nlme")}
  library("nlme")
  
  metaboliteNames = c()
  metadataCatagories = c()
  pvals = c()
  
  #vessel to hold the failures for future reference
  failedLme <<- data.frame(#'<<-' gives global scope
    metaboliteName = character(0), 
    metaboliteRowNum = numeric(0), 
    metadataName = character(0),
    metadataRowNum = numeric(0)
  )
  
  for (i1 in 1:nrow(testingData)){
    #print(paste('i1:',i1))
    for (i2 in 1:nrow(metadata)){
      #print(paste('i2:',i2))
      if (row.names(metadata)[i2] != 'plant number'){
        tryCatch(
          {
            metabo = unlist(testingData[i1,])
            metad = as.factor(unlist(metadata[i2,]))
            pn = unlist(metadata['plant number',])
            myD = data.frame(metabo, metad, pn)
            myLme = lme(metabo ~ metad, 
                        method= "REML", 
                        random = ~1 | pn, 
                        na.action = 'na.omit',
                        data = myD
            )
            myAnova = anova(myLme)
            myPval = myAnova$"p-value"[2]
            metaboliteNames = c(metaboliteNames, row.names(testingData)[i1])
            metadataCatagories = c(metadataCatagories, row.names(metadata)[i2])
            pvals =  c(pvals, myPval)
          },
          error=function(cond) {
            #metaboliteName, metaboliteRowNum, metadataName, metadataRowNum
            myRow = list(row.names(testingData)[i1], i1, row.names(metadata)[i2], i2)
            names(myRow) = c('metaboliteName', 'metaboliteRowNum', 'metadataName', 'metadataRowNum')
            # print(metabo)
            failedLme <<- rbind(failedLme, data.frame(myRow))
            #print('an error is thrown')
            message(cond)
            # Choose a return value in case of error
            # return(NA)
          },
          warning=function(cond) {
            print('a warning is thrown')
            myRow = list(row.names(testingData)[i1], i1, row.names(metadata)[i2], i2)
            names(myRow) = c('metaboliteName', 'metaboliteRowNum', 'metadataName', 'metadataRowNum')
            #print(names(myRow))
            failedLme = rbind(failedLme, data.frame(myRow))
            #print('failedLme warning')
            message(cond)
            # Choose a return value in case of warning
            #return(NULL)
          }
        )
      }
    }
  }
  
  pvalAdj = p.adjust(pvals, method = "BH")
  
  dFrame <- data.frame(pvals,pvalAdj,metaboliteNames,metadataCatagories)
  dFrame <- dFrame[order(dFrame$pvalAdj),]
  
  write.table(dFrame, file=paste(ds, "MetaVs", label, "Lme.txt", sep = ''), row.names=FALSE, sep="\t")
  
  write.table(failedLme, file=paste(ds, "FailedMetaVs", label, "Lme.txt", sep = ''), row.names=FALSE, sep="\t")
  
  pubDf = reshape(dFrame, idvar = 'metaboliteNames', timevar = 'metadataCatagories', direction="wide", sep="", drop = 'pvals')
  write.table(format(pubDf, digits = 3), file=paste(ds, "MetaVs",label,"ReformatShortLme.txt", sep = ''), row.names=FALSE, sep="\t")
  
  pdf(paste(ds, "HistogramMetaVs", label, ".pdf", sep=""))
  hist(pvals, 
       main = paste('Histogram of lme pvalues for ', label),
       xlab = 'Uncorrected pvalues',
       breaks = 50
  )
  dev.off()
  
  pdf(paste(ds, "MetaVs", label, "Lme.pdf", sep=""))
  par(bty = 'l', 
      mar = c(5, 4, 5, 2) + 0.1)
  for( i in 1:nrow(dFrame)){
    aTitle <- paste(  dFrame$metaboliteNames[i], "vs",  dFrame$metadataCatagories[i], "\nAdjusted Pvalue=", dFrame$pvalAdj[i])
    
    plot( as.factor(as.character(metadata[as.character(dFrame$metadataCatagories[i]),])),
          unlist(testingData[as.character(dFrame$metaboliteNames[i]),]),
          main = aTitle,
          xlab = as.character(dFrame$metadataCatagories[i]),
          ylab = as.character(dFrame$metaboliteNames[i]),
          col = as.factor(as.character(metadata[as.character(dFrame$metadataCatagories[i]),]))
    )
  }
  dev.off()
  print('successfull run of lmeMetadata')
}#end function

lmeMetadata(metabolites, 'NoWhole')

# expGroups = unique(unlist(metadata['tissue',]))
# numGroups = length(expGroups)
# iGroups <- vector()
# jGroup <- vector()
# index <- 1
# pValues <- vector()
# metabloliteName <- vector()
# for( m in 1:nrow(mainBioActive) )
# {	
#   for ( i in 1:(numGroups-1) ) 
#   {
#     iData <- mainBioActive[m, which(metadata['tissue', ] == expGroups[i])]
#     
#     for( j in i:numGroups) 
#     {
#       jData <- mainBioActive[m,which(metadata['tissue', ] == expGroups[j])]
#       if( sum(!is.na(iData)) >1 & sum(!is.na(jData)) >1   )
#       {
#         iGroups[index]  <- expGroups[i]
#         jGroup[index] <-expGroups[j]
#         pValues[index] <- t.test(iData, jData)$p.value
#         metabloliteName[index] <- row.names(mainBioActive)[m]
#         index = index + 1
#         #print( paste(index, " " , metabloliteName, "\n"))
#       }
#     }
#   }
# }
# myFrame <-data.frame(metabloliteName,iGroups,jGroup, pValues)
# myFrame$pvalAdj = p.adjust(pValues, method="BH")
# myFrame = myFrame[order(myFrame$pvalAdj),]
# write.table(myFrame, file=paste('lcNeg', 'Mds', lab, 'PairwiseTtest.txt', sep = ''), row.names=FALSE, sep="\t")
# sum(p.adjust(pValues, method="BH") < .1)
# hist(pValues, 
#      breaks = 50,
#      xlab = 'Uncorrected pvalues',
#      main = 'Pairwise t-tests'
# )
