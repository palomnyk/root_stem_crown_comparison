# I want 2 figures for each metadata catagory, 
#   one with plots for each dataset showing all kava lactones with that dataset
#   one with plots for each: Kavain Dihydromethysticin Desmethoxyyangonin
#     for each dataset
# 
# Plan 1: 
#   Read in datasets
#   For metadataRow in metadata rows, 
#     for dataset in datasets:
#       Make boxplot for each dataset showing all kava lactones with that dataset
# Plan 2: 
#   Read in datasets
#   For metadataRow in metadata rows, 
#     for dataset in datasets
#       for bioactive in bioactives
#         Make boxplot for each dataset showing all kava lactones with that dataset

rm(list = ls()) #clear workspace

#set color palette
if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}
palette( brewer.pal(6,"Accent") )

setwd(file.path('~','git','root_stem_crown_comparison'))

#read in meta-data
metadata = read.table(file.path('.','gc','meta_data.csv'), 
                      sep=",", 
                      header=TRUE, 
                      row.names = 1, 
                      check.names = FALSE,
                      stringsAsFactors=FALSE
)

#read in the datasets
gc = read.table(
  file.path('.', 'gc','gc_data_numeric.csv'), 
  sep=",", 
  header=TRUE, 
  row.names = 1, 
  check.names = FALSE
)

lcPos = read.table(
  file.path('.', 'lcPos','lcPos.csv'), 
  sep=",", 
  header=TRUE, 
  row.names = 1, 
  check.names = FALSE
)

lcNeg = read.table(
  file.path('.', 'lcNeg','lcNeg.csv'), 
  sep=",", 
  header=TRUE, 
  row.names = 1, 
  check.names = FALSE
)

datasets = list(gc, lcPos, lcNeg)
datasetNames = c('GC/MS', 'LC/MS+', 'LC/MS-')

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
        mar = c(5, 4, 5, 2) + 0.1)
    plot( as.factor(as.character(metadata[as.character(dFrame$metaName[i]),])),
          unlist(testingData[dFrame$index[i],]),
          main = aTitle,
          xlab = as.character(dFrame$metaName[i]),
          ylab = paste(as.character(dFrame$indexName[i])),
          col = brewer.pal(length(metadata[as.character(dFrame$metaName[i]),]),"Accent")
    )
  }
  dev.off()
  
  print('lmMetadata complete')
}#end of lmMetadata


for (lab in row.names(metadata)){
  print(lab)
  
  
  for(ds in 1:length(datasets)){
    print(datasetNames[ds])
    
    pdf(path.join(getwd(), 'figures', paste(datasetNames[ds], lab, "LinearModFigure.pdf", sep = "")))
    par(par(mfrow=c(1,3)))
    
            myLm = lm( testingData ~ as.factor(as.character(metadata[lab,])))
            metaPval =  c(metaPval, anova(myLm)$"Pr(>F)"[1])

    
  }#end datasets loop
  
}#end row.names(metadata) loop


