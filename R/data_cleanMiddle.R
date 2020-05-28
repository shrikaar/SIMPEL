#pull out the second, or middle, field
data_cleanMiddle <- function(x) sapply (strsplit(x , '[_]' ), `[` , 2)


