#pull out the ending field
data_cleanEnd <- function(x) sapply (strsplit(x , '[_]' ), `[` , 3)
