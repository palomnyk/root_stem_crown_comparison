# Effect size of each dataset
# Effect size is (experimental_mean - control_mean)/stdev
# Current plan:
#   compare gc to lcNeg, gc to lcPos, and lcPos to lcNeg
# 
# pooled stdev = sqrt( ( (N1 - 1)*SD1^2 + (N2 - 1)*SD2^2 ) / (N1 + N2 - 2))

rm(list = ls()) #clear workspace

if (!require("RColorBrewer")) {
  install.packages("RColorBrewer")
  library(RColorBrewer)
}

#set color palette
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
thickness = c(1.75,1.2,0.5)
den = 20
brks = 15
ylm = c(0, 80)

colorz = brewer.pal(length(datasets),"Accent")

rsqrds = list()
metaPval = list()

for (ds in 1:length(datasets)){
  #for row in datasets[ds]
  print('ds')
  print(ds)
  print(colorz[ds])
  dataset = as.data.frame(datasets[ds])
  
  index = c()
  metaName = c()
  metaPval = c()
  indexName = c()
  r_sqrd = c()
  
  for (i1 in 1:nrow(dataset)){
    for (i2 in 1:nrow(metadata)){
      tryCatch(
        { 
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
    pdf(paste("HistogramRsqrdMetadataLm.pdf", sep=""))
    par(lwd=thickness[ds],
        cex.axis = 1.5,
        cex.lab = 1.5)
    hist(r_sqrd,
         main = paste('Histogram lm r squared '),
         xlab = 'Adjusted R squared',
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
         xlab = 'Adjusted R squared',
         breaks = brks,
         col = colorz[ds],
         border = colorz[ds],
         ylim = ylm,
         density = den
    )
  }
}

par(new=TRUE)
legend('topright', 
       legend = datasetNames, 
       col = colorz, 
       pch = 15,
       cex = 1.5
)

dev.off()


