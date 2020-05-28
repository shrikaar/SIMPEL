##function that will create all the dataframes
createDataframe = function(columnnames)
{
  columnnames = columnnames
  print(columnnames)
  lengthColNames = length(columnnames)
  df =  data.frame(matrix(nrow = 0, ncol = lengthColNames))
  colnames(df) = columnnames
  return(df)
}

