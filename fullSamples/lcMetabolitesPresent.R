#Script for making figure showing which metabolites are present in the LC datasets

setwd(file.path('~','git','root_stem_crown_comparison', 'fullSamples'))

lcPos = read.table(
  file.path('.', '..','lcPos','lcPos.csv'), 
  sep=",", 
  header=TRUE, 
  row.names = 1, 
  check.names = FALSE
)

lcNeg = read.table(
  file.path('.', '..','lcNeg','lcNeg.csv'), 
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

for (cm in commonMetabolites){
  lcPos = addEmptyRowIfRnAbsent(cm, lcPos)
  lcNeg = addEmptyRowIfRnAbsent(cm, lcNeg)
}

reorderDf = function(df){
  return(df[ order(row.names(df)),])
  #df[ order(as.numeric(row.names(df))),]
}

lcPos = reorderDf(lcPos)
lcNeg = reorderDf(lcNeg)

write.table(lcPos,
            file = 'lcPosExtended.csv',
            sep = ',',
            row.names=T,
)
write.table(lcNeg,
            file = 'lcNegExtended.csv',
            sep = ',',
            row.names=T)

lcPos[] <- lapply(lcPos, function(x) ifelse(x>0, '+', x))
lcNeg[] <- lapply(lcNeg, function(x) ifelse(x>0, '+', x))
lcPos[is.na(lcPos)] = '-'
lcNeg[is.na(lcNeg)] = '-'

myDf = data.frame(row.names = commonMetabolites)

newColNames = unique(substr(colnames(lcPos),0,3))


combineColumns = function(df){
  myDf = data.frame(row.names = commonMetabolites)
  colIntervals = c(1,2,3)
  for (colName in newColNames){
    #combine columns
    myDf[colName] <- do.call(paste0, df[,colIntervals])
    colIntervals = colIntervals + 3
  }
  return(myDf)
}

lcPos = combineColumns(lcPos)
lcNeg = combineColumns(lcNeg)

write.table(lcPos,
            file = 'lcPosCombinedSymbol.csv',
            sep = ',',
            row.names=T,
)
write.table(lcNeg,
            file = 'lcNegCombinedSymbol.csv',
            sep = ',',
            row.names=T)



