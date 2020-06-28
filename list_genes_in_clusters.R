library("igraph")
library("tidyverse")

# Constants ============================================================

read_from <- "read/updatedTRN.txt" #edge list file to read from
write_to <- "write/data/cluster_gene_list.csv" #file to write the results to

# Helper functions ============================================================

sort_networks <-
    function(list_of_networks) {
        list_of_networks %>%
        map_int(function(network) length(network)) %>%
        magrittr::subtract() %>%
        rank(ties.method = "first") %>%
        magrittr::extract(list_of_networks, .)
    }

# Structure and write the data ============================================================

edgelist <-    
    scan(read_from, skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 
    
trn <-
    edgelist %>%
    graph_from_edgelist(directed = FALSE) 

clusters <- 
    trn %>%
    cluster_fast_greedy(weights=NULL) %>%
    communities() %>%
	sort_networks()

largest_cluster_size <-
	clusters[[1]] %>%
	length()

table <-
	clusters %>%
	map(
		function(vertices) c(vertices, rep(NA, largest_cluster_size - length(vertices))) #add NAs to lists so that each column is the same length
	) %>%
	as_tibble()

table %>% 
	write_csv(write_to, na="", col_names=FALSE)