# New PCA plot

setwd(file.path('~','git','root_stem_crown_comparison'))
outputDir = "3Nov2019"

ds = 'gc'

#read in meta-data
metadata = read.table(file.path('.','data','meta_data_no_whole_noRootVsStem.csv'), 
                      sep=",", 
                      header=TRUE, 
                      row.names = 1, 
                      check.names = FALSE,
                      stringsAsFactors=FALSE)

metabolites = read.table(file.path('.', 'data','gc_data_no_whole.csv'), 
                         sep=",", 
                         header=TRUE, 
                         row.names = 1, 
                         check.names = FALSE)

if(!dir.exists(outputDir)){
  dir.create(outputDir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}
setwd(file.path('~','git','root_stem_crown_comparison', outputDir))

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

pdf('newPCA.pdf')
par(bty = 'l', 
    mar = c(5, 4, 5, 2) + 0.1,
    cex.lab = 1.25
)
plot(myPCA$PC1, myPCA$PC2,
     col='black',
     pch = c(24,21,22)[as.numeric(as.factor(metadata['structural location',]))],
     bg = as.factor(metadata['root vs stem',]),
     cex = 3,
     xlab = 'PCA 1, 74.68%',
     ylab = 'PCA 2, 14.96%',
     main = 'PCA of GC/MS data'
)
plot.new()
legend('topright', 
       legend = c('crown root', 'lateral root', 'stem'),
       pt.bg = as.factor(c('root', 'root', 'stem')),
       pch = c(24,21,22),
       cex = 1
)
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
          col = c(1,1,1,2,2)
    )
  }
  dev.off()
  print('successfull run of lmeMetadata')
}#end function

lmMetadata = function(testingData, label){
  
  index = c()
  metaName = c()
  metaPval = c()
  indexName = c()
  r_sqrd = c()
  
  for (i1 in 1:nrow(testingData)){
    for (i2 in 1:nrow(metadata)){
      tryCatch(
        { 
          myLm = lm( unlist(testingData[i1,]) ~ as.factor(as.character(metadata[i2,])))
          index = c(index, i1)
          metaName = c(metaName, row.names(metadata)[i2])
          indexName = c(indexName, row.names(testingData)[i1])
          metaPval =  c(metaPval, anova(myLm)$"Pr(>F)"[1])
          r_sqrd =  c(r_sqrd, summary(myLm)$r.squared)
        },
        error=function(cond) {
          print('an error is thrown')
          message(cond)
        },
        warning=function(cond) {
          print('a warning is thrown')
          message(cond)
          # Choose a return value in case of warning
          #return(NULL)
        }
      )
    }
  }
  
  pdf(paste(ds, "HistogramRsqrdMetadata", label, "Lm.pdf", sep=""))
  hist(r_sqrd, 
       main = paste(ds, ' Histogram of lm r squared ', label),
       xlab = 'Adjusted R squared',
       breaks = 50
  )
  dev.off()
  
  pdf(paste(ds, "HistogramPvalMetaVs", label, "Lm.pdf", sep=""))
  hist(metaPval, 
       main = paste(ds, ' Histogram of lm pvalues for ', label),
       xlab = 'Uncorrected pvalues',
       breaks = 50
  )
  dev.off()
  
  metaPvalAdj = p.adjust(metaPval, method = "BH")
  dFrame <- data.frame(metaPvalAdj,metaName,index, indexName)
  dFrame <- dFrame [order(dFrame$metaPvalAdj),]
  
  write.table(dFrame, file=paste(ds, "MetaVs",label,"LinearMod.txt", sep = ''), row.names=FALSE, sep="\t")
  
  write.table(format(dFrame, digits = 3), file=paste(ds, "MetaVs",label,"ShortLinearMod.txt", sep = ''), row.names=FALSE, sep="\t")
  
  pubDf = reshape(dFrame, idvar = 'indexName', timevar = 'metaName', direction="wide", sep="", drop = "index")
  write.table(format(pubDf, digits = 3), file=paste(ds, "MetaVs",label,"ReformatShortLinearMod.txt", sep = ''), row.names=FALSE, sep="\t")
  
  #print(pubDf)
  
  pdf(paste(ds,"MetaVs",label,"LinearMod.pdf", sep = ""))
  
  
  for( i in 1:nrow(dFrame)){
    aTitle <- paste(label, ': ', dFrame$indexName[i], " vs ",  dFrame$metaName[i], "\nAdj Pval: ", formatC(dFrame$metaPvalAdj[i]), sep = '')
    
    par(bty = 'l', 
        mar = c(5, 4, 5, 2) + 0.1,
        # cex = 1.3,
        # cex.lab = 1.5,
        cex.axis = 1.45
    )
    plot( as.factor(as.character(metadata[as.character(dFrame$metaName[i]),])),
          unlist(testingData[dFrame$index[i],]),
          main = aTitle,
          xlab = as.character(dFrame$metaName[i]),
          ylab = paste(as.character(dFrame$indexName[i])),
          col = as.factor(as.character(metadata[as.character(dFrame$metaName[i]),]))
    )
    
  }
  dev.off()
  
  print('lmMetadata complete')
}#end of lmMetadata

pairwiseTtest = function(df, metaLvl, lab){
  
  expGroups = unique(unlist(metadata[metaLvl,]))
  numGroups = length(expGroups)
  iGroups <- vector()
  jGroup <- vector()
  index <- 1
  pValues <- vector()
  metabloliteName <- vector()
  for( m in 1:nrow(df) )
  {
    for ( i in 1:(numGroups-1) )
    {
      iData <- df[m, which(metadata[metaLvl, ] == expGroups[i])]
      for( j in i:numGroups)
      {
        jData <- df[m,which(metadata[metaLvl, ] == expGroups[j])]
        if( sum(!is.na(iData)) >1 & sum(!is.na(jData)) >1   )
        {
          iGroups[index]  <- expGroups[i]
          jGroup[index] <-expGroups[j]
          pValues[index] <- t.test(iData, jData)$p.value
          metabloliteName[index] <- row.names(df)[m]
          index = index + 1
        }
      }
    }
  }
  myFrame <-data.frame(metabloliteName,iGroups,jGroup, pValues)
  myFrame$pvalAdj = p.adjust(pValues, method="BH")
  myFrame = myFrame[order(myFrame$pvalAdj),]
  write.table(myFrame, file=paste(ds, lab, metaLvl, 'PairwiseTtest.txt', sep = ''), row.names=FALSE, sep="\t")
  print('PairwiseTtest complete!')
}

lmeMetadata(t(myPCA), 'NoWholePca')
lmeMetadata(metabolites, 'NoWhole')
lmMetadata(metabolites, 'NoWhole')
pairwiseTtest(t(myPCA), 'tissue', 'NoWholePca')
pairwiseTtest(t(myPCA), 'structural location', 'NoWholePca')
# pairwiseTtest(t(myPCA), 'root vs stem', 'NoWholePca')
# pairwiseTtest(metabolites, 'root vs stem', 'NoWholePca')
pairwiseTtest(metabolites, 'tissue', 'NoWholeMetabolites')
pairwiseTtest(metabolites, 'structural location', 'NoWholeMetabolites')
print(summary(group_pca))
