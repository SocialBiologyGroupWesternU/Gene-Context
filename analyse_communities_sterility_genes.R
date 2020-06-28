library("igraph")
library("tidyverse")

# Constants ============================================================

read_from <- "read/updatedTRN.txt" # Edgelist file
write_to <- "write/data/sterility_gene_data_per_cluster.csv" # File to save result data to
sterility_gene_sets_file <- "read/sterility_genes_OGS3_2.csv" # File containing sterility gene sets

#Load and organize base data ============================================================

edgelist <-    
    scan(read_from, skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 

trn <-
    edgelist %>%
    graph_from_edgelist(directed=FALSE) 

clusters <- 
    trn %>%
    cluster_fast_greedy(weights=NULL)

sterility_gene_sets <-
    read_csv(sterility_gene_sets_file) %>%
    gather(key="gene_set", value="gene", na.rm = TRUE)

genes_w_cluster <-
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
	sterility_gene_sets %>%
	inner_join(genes_w_cluster) %>%
	group_by(gene_set,cluster) %>% 
	summarise(gene_count = n())

table %>%	
	write_csv(write_to)
