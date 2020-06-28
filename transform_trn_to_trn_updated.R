library("magrittr")
source("functions.R")

# Constants ============================================================

old_edgelist <- "read/TRN.txt" # AM probe edgelist file
translation_table <- "read/AM_to_OGS3_2.txt" # File that contains AM to OGS 3.2 pairs 
write_to <- "write/data/updatedTRN.txt" # File to write updated edgelist to

# Read and update trn ============================================================

trn <- 
    scan(old_edgelist, skip = 1, what = "character")

am_to_ogs_3_2 <- 
    scan(translation_table,skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 

trn_updated <-
    trn %>%
    translate(am_to_ogs_3_2) %>%
    matrix(ncol=2,byrow=TRUE) %>%
    na.omit() %>%
    unique()

trn_updated %>%
    write.table(write_to, sep = "\t", row.names=FALSE, col.names=FALSE)
