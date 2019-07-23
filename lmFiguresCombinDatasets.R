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


for (lab in row.names(metadata)){
  print(lab)
  print(as.factor(as.character(metadata[lab,])))
  
  pdf(file.path(getwd(), 'figures', paste( lab, "LinearModFigure.pdf", sep = "")))
  par(par(mfrow=c(1,3)))
  
  for(ds in 1:length(datasets)){
    print(datasetNames[ds])
    
    testingData = data.frame(datasets[ds])
    
    catagories = unique(as.character(metadata[lab,]))
    
    testingData[1,]
    
    myRow = as.vector(testingData[1, ])
  
    # It has to be a flattened df -R's nested loops simply wont do it.
    # column 1 == values, column 2 == metadata catagory
    # possible 3rd column == metabolite for bioactives (after subsetting)
    
    lists = list()
    listLevels = levels(as.factor(metadata[lab,]))
    metaFact = as.factor(metadata[lab,])
    metaVect = as.vector(metadata[lab,])
    testDim = dim(testingData)
    dfLength = testDim[1]*testDim[2]
    flatDf = data.frame(numeric(dfLength), character(dfLength))
    
    for (lev in length(listLevels)) {
      print(lev)
      lists[lev] = list(0)
    }
    for(i in 1:dim(testingData)[2]){
      print(match(metaVect[i], listLevels))
      print(length(testingData[,i]))
      index = match(metaVect[i], listLevels)
      append(lists[[index]], testingData[,i])
    }
    newDf = as.data.frame(lists())

    myLm = lm( testingData ~ as.factor(as.character(metadata[lab,])))
    metaPval =  c(metaPval, anova(myLm)$"Pr(>F)"[1])
    
    aTitle <- paste(label, ': ', dFrame$indexName[i], " vs ",  dFrame$metaName[i], "\nAdj Pval: ", formatC(dFrame$metaPvalAdj[i]), sep = '')
    
    par(bty = 'l', 
        mar = c(5, 4, 5, 2) + 0.1)
    plot( as.factor(as.character(metadata[as.character(metadata[lab]),])),
          testingData,
          main = aTitle,
          xlab = as.character(lab),
          ylab = paste(as.character(datasets[ds])),
          #col = brewer.pal(length(metadata[as.character(lab),]),"Accent")
    )
    
  }#end datasets loop
  
}#end row.names(metadata) loop


