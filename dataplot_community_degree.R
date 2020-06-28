library("igraph")
library("tidyverse")
library("scales")

# Constants ============================================================

to_keep <- 7 # Number of clusters to plot in order of size (largest to smallest)
read_from <- "read/updatedTRN.txt" # Edgelist file
to_write <- "write/figure/dataplot/degree_histogram_per_cluster.pdf" # Figure file

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
    group_by(cluster) %>%
    mutate(degree=
        vertices %>%
        induced_subgraph(graph,vids=.) %>%
        degree(.,v=V(.))
    ) %>%
    ungroup() 

#Plot ============================================================

pallete <-
    RColorBrewer::brewer.pal(3,"RdGy")

pdf(file=to_write)

for (i in 1:to_keep)
{
    data <-
        base_data %>%
        filter(
            cluster_rank == i
        )

    plot <- 
        ggplot(data, aes(x=degree)) +
        geom_histogram(aes(fill=type), binwidth=1, alpha=1, position="identity") +
        scale_y_continuous(trans = pseudo_log_trans(base = 10)) + 
        labs(title = paste("Cluster",i, sep = " "),
            x = "Degree",
            y = "Frequency") +
        theme_linedraw() +
        theme(panel.grid = element_blank(), 
            panel.border = element_blank(), 
            axis.line = element_line(colour = "black", size = .2),
            axis.title.x = element_text(size=20),
            axis.title.y = element_text(size=20),
            legend.title = element_text(size=20),
            title = element_text(size=22),
            legend.position = "none") +
        scale_fill_manual(values=pallete[c(3,1)]) +
        xlim(0, 160) 

    plot %>%
        print()
}

dev.off()