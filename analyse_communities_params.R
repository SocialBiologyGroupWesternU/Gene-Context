library("igraph")
library("tidyverse")

# Constants ============================================================

vertex_count_threshold <- 10 # Any clusters with a vertex count less than this will be dropped
read_from <- "read/updatedTRN.txt" # Edgelist file
write_to <- "write/data/assorted_parameters_per_cluster.csv" # File to save result data to

# Helper functions ============================================================

apply_function_per_network <-
	function(fun) {
		networks_to_analyse %>%
		map(function(network) fun(network)) %>%
        unlist()
	}

sort_networks <-
    function(list_of_networks) {
        list_of_networks %>%
        map_int(function(network) gorder(network)) %>%
        magrittr::subtract() %>%
        rank(ties.method = "first") %>%
        magrittr::extract(list_of_networks, .)
    }

#Load and organize base data ============================================================

edgelist <-    
    scan(read_from, skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 
    
trn <-
    edgelist %>%
    graph_from_edgelist(directed = FALSE) 

networks_to_analyse <- 
    trn %>%
    cluster_fast_greedy(weights=NULL) %>%
    communities() %>%
    map(
        function(cluster) induced_subgraph(trn, cluster)
    ) %>%
    keep(
        function(network) gorder(network) > vertex_count_threshold #discard networks that dont meet the degree threshold
    ) %>%
    sort_networks() %>%
    c(list(trn)) #add trn to the list so it is included in the analysis

# Build data table from base data ============================================================

table <-
    tibble(number_of_nodes=
        apply_function_per_network(igraph::gorder)
    ) %>%
    mutate(transcription_factors=
        apply_function_per_network(
            function(network){
                network %>% 
                V() %>% 
                as_ids() %>% 
                intersect(edgelist[,1]) %>% 
                length()
            }
        )
    ) %>%
    mutate(target_genes=
        apply_function_per_network(
            function(network){
                network %>% 
                V() %>% 
                as_ids() %>% 
                intersect(edgelist[,2]) %>% 
                length()
            }
        )
    ) %>%
    mutate(number_of_links=
        apply_function_per_network(igraph::gsize)
    ) %>%
    mutate(average_degree=
        apply_function_per_network(
            function(network) {
                network %>%
				degree() %>%
                mean()
			}
        )
    ) %>%
    mutate(k_max=
    	apply_function_per_network(
            function(network) {
                network %>%
				degree() %>%
                max()
			}
        )
    ) %>%
    mutate(hub_gene=
        apply_function_per_network(
            function(network) {
                vertices <- network %>% V() %>% as_ids()
                network %>%
				degree(v=vertices) %>%
                order(decreasing=TRUE) %>%
                first() %>%
                nth(vertices,n=.)
			}
        )
    ) %>%
    mutate(mean_betweenness=
        apply_function_per_network(
            function(network) {
                network %>%
				betweenness() %>%
                mean()
			}
        )
    ) %>%
    mutate(mean_internode_distance=
    	apply_function_per_network(igraph::mean_distance)
    ) %>%
    mutate(assortativity=
    	apply_function_per_network(igraph::assortativity_degree)
    ) %>%
    mutate(k1_count=
        apply_function_per_network(
            function(network) {
                network %>%
				degree() %>%
                keep(function(deg) deg == 1) %>%
                length()
			}
        )
    )%>%
    mutate(alpha=
        apply_function_per_network(
            function(network) {
                network %>%
				degree() %>%
                fit_power_law() %>%
                magrittr::use_series(alpha)
			}
        )
    ) %>%
    mutate(average_degree_tf=
        apply_function_per_network(
            function(network) {
                network %>% 
                V() %>% 
                as_ids() %>% 
                intersect(edgelist[,1]) %>%
                degree(network, v=.) %>%
                mean()
			}
        )
    )%>%
    mutate(average_degree_tg=
        apply_function_per_network(
            function(network) {
                network %>% 
                V() %>% 
                as_ids() %>% 
                intersect(edgelist[,2]) %>%
                degree(network, v=.) %>%
                mean()
			}
        )
    )%>%
    mutate(motifs_3=
        apply_function_per_network(
            function(network) count_motifs(network, size=3)
        )
    )%>%
    mutate(motifs_4=
        apply_function_per_network(
            function(network) count_motifs(network, size=4)
        )
    )

table %>%
    write_csv(write_to)