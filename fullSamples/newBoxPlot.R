# New PCA plot

setwd(file.path('~','git','root_stem_crown_comparison'))

#read in meta-data
metadata = read.table(file.path('.','data','meta_data_no_whole.csv'), 
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


lmeMetadata = function(testingData){
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
  
  write.table(dFrame, file="NewPlotsMetaVsLme.txt", row.names=FALSE, sep="\t")
  
  pdf("NewPlotsMetaVsLme.pdf")
  par(bty = 'l', 
    mar = c(5, 4, 5, 2) + 0.1
    )
  for( i in 1:nrow(dFrame)){
    if (dFrame$metadataCatagories[i] == 'tissue' & dFrame$metaboliteNames[i] %in% c('Kavain', 'Dihydromethysticin', 'Desmethoxyyangonin')){
      aTitle <- paste(  dFrame$metaboliteNames[i], "vs",  dFrame$metadataCatagories[i], "\nAdjusted Pvalue=", dFrame$pvalAdj[i])
      
      plot( as.factor(as.character(metadata[as.character(dFrame$metadataCatagories[i]),])),
            unlist(testingData[as.character(dFrame$metaboliteNames[i]),]),
            main = aTitle,
            xlab = as.character(dFrame$metadataCatagories[i]),
            ylab = as.character(dFrame$metaboliteNames[i]),
            col = as.factor(c('root', 'root', 'root', 'stem', 'stem'))
      )
      legend('topright', 
             legend = c('root', 'stem'),
             fill = as.factor(c('root', 'stem')),
             cex = 1
      )
    }
    
  }
  dev.off()
  print('successfull run of lmeMetadata')
}#end function

lmeMetadata(metabolites)