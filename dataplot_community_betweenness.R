library("igraph")
library("tidyverse")

# Constants ============================================================

to_keep <- 7 # Number of clusters to plot in order of size (largest to smallest)
read_from <- "read/updatedTRN.txt" # Edgelist file
to_write <- "write/figure/dataplot/betweenness_boxplot_all_clusters.pdf" # Figure file

#Load Data ============================================================

edgelist <-    
    scan(read_from, skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 

graph <-    
    edgelist %>%
    graph_from_edgelist(directed = FALSE) 

clusters_object <- 
    graph %>%
    cluster_fast_greedy(weights=NULL)	

#Build table to plot ============================================================

data <-     
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
    group_by(cluster) %>%
    mutate(betweenness=
        vertices %>%
        induced_subgraph(graph,vids=.) %>%
        betweenness(.,v=V(.), normalized=TRUE)
        ) %>%
    ungroup() %>%
    filter(
        cluster_rank <= to_keep
    )

#Plot ============================================================

pallete <-
    RColorBrewer::brewer.pal(3,"RdGy")

plot <- 
    ggplot(data, aes(x=cluster_rank %>% factor(), y=betweenness, fill=type)) +
    geom_boxplot() +
    labs(x = "Cluster",
        y = "Betweeness",
        fill = "Gene Type") +
    theme_linedraw() +
    theme(panel.grid = element_blank()) +
    scale_fill_manual(values=pallete[c(3,1)]) +
    coord_cartesian(ylim = c(0,0.4))

ggsave(to_write)

