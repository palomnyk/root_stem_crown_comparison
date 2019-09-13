#Script for making figure showing which metabolites are present in the LC datasets

setwd(file.path('~','git','root_stem_crown_comparison'))

lcPos = read.table(
  file.path('.', 'data','lcPosNoWhole.csv'), 
  sep=",", 
  header=TRUE, 
  row.names = 1, 
  check.names = FALSE
)

lcNeg = read.table(
  file.path('.', 'data','lcNegNoWhole.csv'), 
  sep=",", 
  header=TRUE, 
  row.names = 1, 
  check.names = FALSE
)

commonMetabolites = sort(unique(c(row.names(lcPos), row.names(lcNeg))))

addEmptyRowIfRnAbsent = function(rowName, df){
  if (! rowName %in% row.names(df)){
    emptyVec = rep(NA, ncol(df))
    names(emptyVec) = colnames(df)
    empty = data.frame(emptyVec)
    df = rbind(df, t(empty))
    rownames(df)[nrow(df)] = cm
    return(df)
  }
  else{
    return(df)
  }
}

lcPos[] <- lapply(lcPos, function(x) ifelse(x>0, '+', x))
lcNeg[] <- lapply(lcNeg, function(x) ifelse(x>0, '+', x))

for (cm in commonMetabolites){
  lcPos = addEmptyRowIfRnAbsent(cm, lcPos)
  lcNeg = addEmptyRowIfRnAbsent(cm, lcNeg)
}

lcPos[is.na(lcPos)] = '-'
lcNeg[is.na(lcNeg)] = '-'

myDf = data.frame(row.names = commonMetabolites)

newColNames = unique(substr(colnames(lcPos),0,3))

colIntervals = c(1,2,3)
for (colName in newColNames){
  #combine columns
  myDf[colName] <- do.call(paste0, lcPos[,colIntervals])

  colIntervals = colIntervals + 3
  print(colIntervals)
}

print(newColNames)

write.table(lcPos,
            file = 'lcPosExtended.csv',
            sep = ',',
            row.names=T,
            )
write.table(lcNeg,
            file = 'lcNegExtended.csv',
            sep = ',',
            row.names=T)


