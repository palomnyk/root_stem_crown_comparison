# Script for creating tiffs for giga science

setwd(file.path('~','git','root_stem_crown_comparison'))
ds = 'lcPos'
outputDir = file.path('18Aug2020', ds)

#read in metadata
metadata = read.table(file.path('.','data','meta_data_no_whole_noRootVsStem.csv'), 
                      sep=",", 
                      header=TRUE, 
                      row.names = 1, 
                      check.names = FALSE,
                      stringsAsFactors=FALSE)

old_metadata = read.table(file.path('.','data','meta_data_no_whole.csv'), 
                          sep=",", 
                          header=TRUE, 
                          row.names = 1, 
                          check.names = FALSE,
                          stringsAsFactors=FALSE)

metabolites = read.table(file.path('.', 'data', paste0(ds, '_data_no_whole.csv')), 
                         sep=",", 
                         header=TRUE, 
                         row.names = 1, 
                         check.names = FALSE)

if(!dir.exists(outputDir)){
  dir.create(outputDir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}
if(!dir.exists(file.path(outputDir, "PCA"))){
  dir.create(file.path(outputDir, "PCA"), showWarnings = TRUE, recursive = FALSE, mode = "0777")
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

tiff(file.path( 'PCA', paste0(ds, 'PCA.tiff')),
     width = 1200, height = 1200, units = "px", pointsize = 12,
     compression = "none",
)
par(bty = 'l', 
    mar = c(5, 4, 5, 2) + 0.1,
    cex.lab = 1.25
)
plot(myPCA$PC1, myPCA$PC2,
     col='black',
     pch = c(24,21,22)[as.numeric(as.factor(metadata['structural location',]))],
     bg = as.factor(old_metadata['root vs stem',]),
     cex = 3,
     #this part needs to be customized for each figure
     xlab = 'PCA1: 36%',
     ylab = 'PCA2: 21%',
     main = 'PCA of LC-MS Positive data'
)
dev.off()
tiff(file.path( 'PCA', paste0(ds, 'PCALegend.tiff')),
     width = 400, height = 400, units = "px", pointsize = 12,
     compression = "none",
)
plot.new()
legend('topright', 
       legend = c('crown root', 'lateral root', 'stem'),
       pt.bg = as.factor(c('root', 'root', 'stem')),
       pch = c(24,21,22),
       cex = 1
)
dev.off()

lmMetadata = function(testingData, label){#this is the example
  
  index = c()
  metaName = c()
  metaPval = c()
  indexName = c()
  r_sqrd = c()
  
  for (i1 in 1:nrow(testingData)){
    for (i2 in 1:nrow(metadata)){
      tryCatch(
        { 
          myLm = lm( unlist(testingData[i1,]) ~ as.factor(as.character(metadata[i2,])))#factor are an easy way to subset your data in R
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
  
  metaPvalAdj = p.adjust(metaPval, method = "BH")#adjusting pvalues - you could skip this
  dFrame <- data.frame(metaPvalAdj,metaName,index, indexName)
  dFrame <- dFrame [order(dFrame$metaPvalAdj),]
  
  for( i in 1:nrow(dFrame)){
    
    aTitle <- paste(label, ': ', dFrame$indexName[i], " vs ",  dFrame$metaName[i], "\nAdj Pval: ", formatC(dFrame$metaPvalAdj[i]), sep = '')
    
    tiff(paste(ds,dFrame$indexName[i], "vs",  dFrame$metaName[i],label,"LinearMod.tiff", sep = "_"),
         width = 1200, height = 1200, units = "px", pointsize = 12,
         compression = "none",
    )
    par(bty = 'l', 
        mar = c(5, 4, 5, 2) + 0.1,
        cex.axis = 1.45
    )
    plot( as.factor(as.character(metadata[as.character(dFrame$metaName[i]),])),
          unlist(testingData[dFrame$index[i],]),
          main = aTitle,
          xlab = as.character(dFrame$metaName[i]),
          ylab = paste(as.character(dFrame$indexName[i])),
          col = as.factor(as.character(metadata[as.character(dFrame$metaName[i]),]))
    )
    dev.off()
  }
  
  print('lmMetadata complete')
}#end of lmMetadata
lmMetadata(metabolites, 'NoWhole')