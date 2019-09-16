# Effect size of each dataset
# Effect size is (experimental_mean - control_mean)/stdev
# Current plan:
#   compare gc to lcNeg, gc to lcPos, and lcPos to lcNeg
# 
# pooled stdev = sqrt( ( (N1 - 1)*SD1^2 + (N2 - 1)*SD2^2 ) / (N1 + N2 - 2))

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

#create PCA
gc_pca = data.frame(prcomp(t(na.omit(gc)), 
                   center = TRUE,
                   scale = TRUE)$x
)

lcPos = read.table(
  file.path('.', 'lcPos','lcPos.csv'), 
  sep=",", 
  header=TRUE, 
  row.names = 1, 
  check.names = FALSE
)

#create PCA
lcPos_pca = data.frame(prcomp(t(na.omit(lcPos)), 
                           center = TRUE,
                           scale = TRUE)$x
)

lcNeg = read.table(
  file.path('.', 'lcNeg','lcNeg.csv'), 
  sep=",", 
  header=TRUE, 
  row.names = 1, 
  check.names = FALSE
)

#create PCA
lcNeg_pca = data.frame(prcomp(t(na.omit(lcNeg)), 
                           center = TRUE,
                           scale = TRUE)$x
)

datasets = list(gc_pca, lcPos_pca, lcNeg_pca)
datasetNames = c('GC/MS', 'LC/MS+', 'LC/MS-')
thickness = c(1.5,1,0.5)
den = 20
brks = 'Scott'
ylm = c(0, 120)

colorz = brewer.pal(length(datasets),"Accent")

rsqrds = data.frame(gcms = numeric(), 
                    lcmsP = numeric(), 
                    lcmsN = numeric()
)
metaPval = list()

for (ds in 1:length(datasets)){
  #for row in datasets[ds]
  print('ds')
  dataset = as.data.frame(datasets[ds])
  print(dim(dataset))
  
  index = c()
  metaName = c()
  metaPval = c()
  indexName = c()
  r_sqrd = c()
  
  for (i1 in 1:nrow(dataset)){
    for (i2 in 1:nrow(metadata)){
      tryCatch(
        { 
          print(row.names(dataset)[i1])
          print(length(dataset[i1,]))
          print(metadata[i2,])
          p
          myLm = lm( unlist(dataset[i1,]) ~ as.factor(as.character(metadata[i2,])))
          index = c(index, i1)
          metaName = c(metaName, row.names(metadata)[i2])
          indexName = c(indexName, row.names(dataset)[i1])
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
  
  if (ds == 1){
    pdf(paste("HistogramRsqrdMetadataPCALm.pdf", sep=""))
    par(lwd=thickness[ds],
        cex.axis = 1.5,
        cex.lab = 1.5)
    hist(r_sqrd,
         main = paste('Histogram lm r squared '),
         xlab = 'R squared',
         breaks = brks,
         col = colorz[ds],
         border = colorz[ds],
         ylim= ylm,
         density = den
    )
    count <- table(r_sqrd)  
  }else{
    par(new=TRUE,
        lwd=thickness[ds])
    hist(r_sqrd,
         main = paste('Histogram lm r squared '),
         xlab = 'R squared',
         breaks = brks,
         col = colorz[ds],
         border = colorz[ds],
         ylim = ylm,
         density = den
    )
    par(new=TRUE)
    legend('topright', 
           legend = datasetNames, 
           col = colorz, 
           pch = 15,
           cex = 1.5
    )
  }
}

dev.off()
