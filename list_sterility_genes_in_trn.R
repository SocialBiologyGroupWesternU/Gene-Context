library("igraph")
library("tidyverse")

# Constants ============================================================

trn_file <- "read/updatedTRN.txt" # Edgelist file
sterility_gene_sets_file <- "read/sterility_genes_OGS3_2.csv" # Sterility genes file
write_to <- "write/data/sterility_genes_trn.csv" # File to save result data to

# Load, organize and write the data ============================================================

edgelist <-    
    scan(trn_file, skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 
    
trn <-
    edgelist %>%
    graph_from_edgelist(directed = FALSE) 

clusters <- 
    trn %>%
    cluster_fast_greedy(weights=NULL)

sterility_gene_sets <-
    read_csv(sterility_gene_sets_file) %>%
    gather(key="gene_set", value="gene", na.rm = TRUE)

trn_genes_w_cluster <-
    tibble(
        gene=
            trn %>%
            V() %>%
            as_ids(),
        cluster=
            clusters %>% 
            membership() %>%
            magrittr::extract(
                clusters %>%
                sizes() %>%
                magrittr::subtract() %>%
                rank(ties.method = "first"),
                .
            )
    )

table <-
    inner_join(sterility_gene_sets,trn_genes_w_cluster) %>%
    arrange(gene_set,cluster)

table %>%
    write_csv(write_to)