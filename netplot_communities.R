library("igraph")

#Constants ============================================================

to_keep <- 9 # Number of clusters to plot in order of size (largest to smallest)
read_from <- "read/updatedTRN.txt" # Edgelist file
write_to <- "write/figure/netplot/clusters_bipartite_layout.pdf" # Figure file

#Load and organize the data ============================================================

edgelist <-    
    scan("read/updatedTRN.txt", skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 

graph <-    
    edgelist %>%
    graph_from_edgelist(directed=FALSE) 

clusters_object <- 
	graph %>%
	cluster_fast_greedy(weights=NULL)

ordered_clusters <-
	clusters_object %>%
	communities() %>% 
	magrittr::extract(
		clusters_object %>%
		sizes() %>%
		order(decreasing=TRUE)
	)

TF_vector <- 
	edgelist[,1]

GENE_vector <- 
	edgelist[,2]

#Plot ============================================================

node_pallete <-
	RColorBrewer::brewer.pal(3,"RdGy")

pdf(file="write/figure/netplot/clusters_kk_layout.pdf") 	

for (i in 1:size)
{
	cluster <- 
		graph %>%
		induced_subgraph(
			ordered_clusters[i] %>%
			unlist()
		)

	vertices <- 
		cluster %>%
		V()

	edges <- 
		cluster %>%
		E()	

	hub_node <- 
		vertices %>%
		magrittr::extract(
			cluster %>%
			degree() %>%
			which.max()
		)

	hub_node_adjacent_edges <-
		edges %>%
		magrittr::extract(
			hub_node %>%
			.inc()
		)
	
	remaining_edges <- 
		edges %>%
		setdiff(
			hub_node_adjacent_edges
		)

	tf_nodes <-
		vertices %>%
		as_ids() %>% 
		magrittr::is_in(
			TF_vector
		)
	
	gene_nodes <- 
		vertices %>%
		as_ids() %>% 
		magrittr::is_in(
			GENE_vector
		)

	degrees <-
		cluster %>%
		degree(vertices) %>%
		magrittr::divide_by(
			max(.) %>%
			magrittr::divide_by(7)
		) %>%
		replace(
			magrittr::is_less_than(.,1.5),
			1.5
		)

	edge_attr(cluster, "color", hub_node_adjacent_edges) <- "Grey35"
	edge_attr(cluster,"color", remaining_edges) <- "Grey85"
	vertex_attr(cluster,"color", tf_nodes) <- node_pallete[1]
	vertex_attr(cluster,"color", gene_nodes) <- node_pallete[3]
	V(cluster)$size <- degrees
	vertex_attr(cluster,"label", vertices) <- ""
	graph_attr(cluster, "layout") <- layout_with_kk

	cluster %>%
		plot(main=paste("Cluster",i, sep=" ")) %>%
		print()
}

dev.off()

