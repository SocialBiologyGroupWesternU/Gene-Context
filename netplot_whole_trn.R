library("igraph")
library("tidyverse")
library("RColorBrewer")

# Constants ============================================================

read_from <- "read/updatedTRN.txt" # Edgelist file
write_to <- "write/figure/netplot/whole_trn_kk_layout.pdf" # Figure file

# Helper functions ============================================================	

alpha <- 
	function(edge){
		incident_vertices <- ends(trn, edge)
		vertex1_deg <- trn %>% degree(v=incident_vertices[1,1])
		vertex2_deg <- trn %>% degree(v=incident_vertices[1,2])
		ifelse(vertex1_deg >= vertex2_deg, vertex1_deg, vertex2_deg) %>% 
		magrittr::multiply_by(255/max_degree) %>% 
		trunc() %>%
		as.hexmode() %>%
		format(width=2, upper.case = TRUE) %>%
		as.character()
	}

colour <- 
	function(edge){
		for (i in 1:length(clusters)){
			if (edge %in% as_ids(E(clusters[[i]])))
				return(colours[i])
		}
		return(NULL)
	}

#Load and organize the data ============================================================

edgelist <-    
    scan(read_from, skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 
    
trn <-
    edgelist %>%
    graph_from_edgelist(directed = FALSE) 

trn_clustered <-
	trn %>%
	cluster_fast_greedy(weights=NULL)

clusters <-
	trn_clustered %>%
	communities() %>%
	map(function(cluster) induced_subgraph(trn, cluster))

colours <- brewer.pal(length(clusters),"Set3")

max_degree <-
	trn %>%
	degree() %>%
	max()

#Add visual attributes to trn then plot ============================================================

trn_with_attr <-
	trn %>%
	set_vertex_attr("label", value="") %>%
	set_vertex_attr("size", value=
		trn %>%
		degree(normalized=TRUE) %>%
		magrittr::multiply_by(30)
	) %>%
	set_vertex_attr("color", value=
		trn_clustered %>%
		membership() %>%
		as.integer() %>%
		map_chr(function(cluster_num) colours[cluster_num])
	) %>%
	set_vertex_attr("frame.color", value=NA) %>%
	set_edge_attr("color", value=
		trn %>%
		E() %>%
		as_ids() %>%
		map_chr(function(edge){
			edge_colour <- colour(edge)
			if(!is.null(edge_colour))
				paste(edge_colour,alpha(edge), sep="")
			else
				NA
			}	
		)
	) %>%
	set_edge_attr("width", value=.5) %>%
	set_graph_attr("layout", layout_with_kk)

pdf(file=write_to) 	
trn_with_attr %>% plot() %>% print()
dev.off()
