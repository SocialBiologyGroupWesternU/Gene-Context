library("igraph")
library("tidyverse")
source("functions.R")

#Load and Organize Data ============================================================

edgelist <-    
    scan("read/updatedTRN.txt", skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 

graph <-    
    edgelist %>%
    graph_from_edgelist(directed = FALSE) 

clusters_object <- 
    graph %>%
    cluster_fast_greedy(weights=NULL)

# Constants ============================================================

to_keep <- 7

pallete <-
    RColorBrewer::brewer.pal(3,"RdGy")
    
#Build Base Data Frame ============================================================

base_data <-     
    graph %>%
    V() %>%
    as_ids() %>%
    enframe(name=NULL, value="vertices") %>%
    mutate(cluster=
        clusters_object %>%
        membership() 
    ) %>%
    mutate(cluster_rank=
        clusters_object %>%
        sizes() %>%
        magrittr::subtract() %>%
        rank(ties.method = "first") %>%
        magrittr::extract(
            cluster
        ) 
    ) %>%
    mutate(type=
        if_else(
            vertices %>% magrittr::is_in(edgelist[,1]),
            "Transcription Factor",
            "Target"
        )
    ) %>%
    filter(
        cluster_rank <= to_keep  
    )

am_to_ogs_3_2 <- 
    scan("read/AM_to_OGS3_2.txt",skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 

#Helper Functions ============================================================

read_gene_set <-
	function(gene_set_location)
	{
		gene_set_location %>%
		scan(skip = 1, what = "character") %>%
    	translate(am_to_ogs_3_2) %>%
    	na.omit()
	}

extend_base_data <-
    function(gene_set)
    {
        base_data %>% 
        mutate(sterility_gene=
            if_else(
                vertices %>% magrittr::is_in(gene_set),
                "Sterility Genes",
                "Remaining Genes"
            )
        )
    }

plot_extended_data <-
    function(extended_data,title)
    {
        pallete <-
            RColorBrewer::brewer.pal(3,"RdGy")
        ggplot(extended_data, aes(x=cluster_rank %>% factor(),y=sterility_gene,colour=type)) +
        geom_jitter() +
        geom_vline(xintercept = c(1:6) + .5, size=.1) +
        labs(title = title,
            x = "Clusters",
            y = NULL,
            colour = "Gene Type") +
        scale_colour_manual(values=pallete[c(3,1)])+
        theme_linedraw() + 
        theme(
            panel.grid.major = element_blank() 
        )
    }

#Gene Set Plotting ============================================================

pdf(file="write/figure/dataplot/sterility_genes_per_cluster.pdf") 

"read/Cardoen_2011_Sterility_Genes.txt" %>%
read_gene_set() %>%
extend_base_data() %>%
plot_extended_data("Cardoen 2011 Gene Set")

print(plot)

"read/Grozinger_2003_Sterility_Genes.txt" %>%
read_gene_set() %>%
extend_base_data() %>%
plot_extended_data("Grozinger 2003 Gene Set")

print(plot)

"read/Grozinger_2007_Sterility_Genes.txt" %>%
read_gene_set() %>%
extend_base_data() %>%
plot_extended_data("Grozinger 2007 Get Set")

print(plot)

"read/Mullen_2014_Sterility_Genes.txt" %>%
read_gene_set() %>%
extend_base_data() %>%
plot_extended_data("Mullen 2014 Gene Set")

print(plot)

"read/Galbraith_2016_Sterility_Genes.txt" %>%
scan(skip = 1, what = "character") %>%
extend_base_data() %>%
plot_extended_data("Galbraith 2016 Gene Set")

print(plot)

dev.off()