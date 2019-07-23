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
flatDatasets = list()

interest = c('tissue', 'structural location', 'root vs stem')


pvals = c()
catagory = c()
datasetNum = c()

pdf(file.path(getwd(), 'figures', paste("fullDsLm.pdf", sep = "")))
#par(mfrow=c(3,3)) #uncomment when using interst in first loop

for (lab in row.names(metadata)){
  print(lab)
  
  for(ds in 1:length(datasets)){
    print(datasetNames[ds])
    
    testingData = log2(data.frame(datasets[ds]) + 1)
    
    #make new flat df
    origDim = dim(testingData)
    totalCells = origDim[1]*origDim[2]
    metabol = character(totalCells)
    auc = numeric(totalCells)
    metaDat = character(totalCells)
    appendCounter = 1
    appendEnd = origDim[1]
    for(i in 1:origDim[2]){
      metabol[appendCounter:appendEnd] = row.names(testingData)
      metaDat[appendCounter:appendEnd] = rep(metadata[lab,i], origDim[1])
      auc[appendCounter:appendEnd] = testingData[,i]
      appendCounter = appendCounter + origDim[1]
      appendEnd = appendEnd + origDim[1]
    }
    flatDf = data.frame(metabol, auc, metaDat)
    flatDatasets[ds] = flatDf;
    
    myLm = lm( flatDf$auc ~ as.factor(flatDf$metaDat))
    metaPval =  anova(myLm)$"Pr(>F)"[1]
    print(metaPval)

    aTitle <- paste(lab, datasetNames[ds], "\nPval:", formatC(metaPval), sep = ' ')

    par(bty = 'l',
        mar = c(1, 1, 1, 1) + 1)
    plot( flatDf$auc ~ as.factor(flatDf$metaDat),
          main = aTitle,
          xlab = as.character(lab),
          ylab = ''
          #unique(as.character(metadata[lab,]))
          #col = brewer.pal(length(metadata[as.character(lab),]),"Accent")
    )
    stripchart(flatDf$auc ~ as.factor(flatDf$metaDat),
               method='jitter',
               jitter=.2,
               vertical=TRUE,
               add=T,
               pch=16,
               cex = 0.6,
               col = 'grey'
    )

  }#end datasets loop

  
}#end row.names(metadata) loop

dev.off()

