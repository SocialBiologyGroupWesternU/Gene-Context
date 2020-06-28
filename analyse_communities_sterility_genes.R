library("igraph")
library("tidyverse")
source("functions.R")

# Constants ============================================================

vertex_count_threshold <- 10 # Any clusters with a vertex count less than this will be dropped
read_from <- "read/updatedTRN.txt" # Edgelist file
write_to <- "write/data/sterility_gene_data_per_cluster.csv" # File to save result data to

# Helper functions ============================================================

read_gene_set <-
	function(gene_set_location) {
		gene_set_location %>%
		scan(skip = 1, what = "character") %>%
    	translate(am_to_ogs_3_2) %>%
    	na.omit()
	}

intersection_per_cluster <-
	function(to_intersect) {
		clusters %>%
		map_int(
			function(cluster) {
				cluster %>%
				intersect(to_intersect) %>%
				length()
			}
		)
	}

#Load and organize base data ============================================================

edgelist <-    
    scan(read_from, skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 

graph <-    
    edgelist %>%
    graph_from_edgelist(directed = FALSE) 

clusters <- 
    graph %>%
    cluster_fast_greedy(weights=NULL) %>%
	communities() %>%
	keep(function(cluster) length(cluster) > vertex_count_threshold) #Discard clusters with a vertex count less than the threshold

am_to_ogs_3_2 <- 
    scan("read/AM_to_OGS3_2.txt",skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 

#Load sterility gene data ============================================================

cardoen_2011_genes <- 
    read_gene_set("read/Cardoen_2011_Sterility_Genes.txt")

grozinger_2003_genes <- 
    read_gene_set("read/Grozinger_2003_Sterility_Genes.txt")

grozinger_2007_genes <- 
    read_gene_set("read/Grozinger_2007_Sterility_Genes.txt")

mullen_2014_genes <- 
    read_gene_set("read/Mullen_2014_Sterility_Genes.txt")

# Build data table ============================================================

table <-
	tibble(number_of_genes=
		clusters %>%
		map_int(length)
	) %>%
	mutate(target_gene_count=
		intersection_per_cluster(edgelist[,2])
	) %>%
	mutate(tf_count=
		intersection_per_cluster(edgelist[,1])
	) %>%
	mutate(cardoen_2011_count=
		intersection_per_cluster(cardoen_2011_genes)
	) %>%
	mutate(grozinger_2003_count=
		intersection_per_cluster(grozinger_2003_genes)
	) %>%
	mutate(grozinger_2007_count=
		intersection_per_cluster(grozinger_2007_genes)
	) %>%
	mutate(mullen_2014_count=
		intersection_per_cluster(mullen_2014_genes)
	) %>%
	arrange(
		desc(number_of_genes)
	)

table %>%	
	write_csv(write_to)
