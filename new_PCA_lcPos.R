# New PCA plot

setwd(file.path('~','git','root_stem_crown_comparison'))

#read in meta-data
metadata = read.table(file.path('.','data','meta_data_no_whole.csv'), 
                      sep=",", 
                      header=TRUE, 
                      row.names = 1, 
                      check.names = FALSE,
                      stringsAsFactors=FALSE)

metabolites = read.table(file.path('.', 'data','lcPosNoWhole.csv'), 
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

pdf('newPCA_lcPos.pdf')
par(bty = 'l', 
    mar = c(5, 4, 5, 2) + 0.1,
    cex.lab = 1.25
)
plot(myPCA$PC1, myPCA$PC2,
     col='black',
     pch = c(24,21,22)[as.numeric(as.factor(metadata['structural location',]))],
     bg = as.factor(metadata['root vs stem',]),
     cex = 3,
     xlab = 'PCA 1, 36.15%',
     ylab = 'PCA 2, 20.82%',
     main = 'PCA of LC/MS Positive data'
)
# legend('topright', 
#        legend = c('crown root', 'lateral root', 'stem'),
#        pt.bg = as.factor(c('root', 'root', 'stem')),
#        pch = c(24,21,22),
#        cex = 1
# )
dev.off()