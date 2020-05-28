#pull out the first field
data_clean <- function(x) sapply (strsplit(x , '[_]' ), `[` , 1)
