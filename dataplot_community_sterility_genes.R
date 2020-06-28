library("igraph")
library("tidyverse")

# Constants ============================================================

to_keep <- 7 # Number of clusters to plot in order of size (largest to smallest)
trn_file <- "read/updatedTRN.txt" #edge list file to read from
sterility_gene_sets_file <- "read/sterility_genes_OGS3_2.csv"
write_to <- "write/figure/dataplot/sterility_genes_per_cluster.pdf" #file to write figures to

# Collect and plot the data ============================================================

edgelist <-    
    scan(trn_file, skip = 1, what = "character") %>%
    matrix(ncol=2,byrow=TRUE) 

trn <-
    edgelist %>%
    graph_from_edgelist(directed=FALSE) 

clusters <- 
    trn %>%
    cluster_fast_greedy(weights=NULL)

gene_w_gene_set <-
    read_csv(sterility_gene_sets_file) %>%
    gather(key="sterility_gene_set", value="gene", na.rm = TRUE)

genes_w_type <-
    tibble(gene=edgelist[,1]) %>%
    mutate(type=rep("Transcription Factor", length(edgelist[,1]))) %>%
    dplyr::union(
        tibble(gene=edgelist[,2]) %>%
        mutate(type=rep("Target Gene", length(edgelist[,2])))
    )

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
    ) %>%
    filter(
        cluster <= to_keep
    )

genes_w_type_w_cluster <-
    genes_w_cluster %>%
    inner_join(genes_w_type) 

gene_sets <-
    gene_w_gene_set %>% 
    pull(sterility_gene_set) %>% 
    unique()

#Plot ============================================================

pallete <- RColorBrewer::brewer.pal(3,"RdGy")

pdf(file=write_to) 	

for (gene_set in gene_sets) 
{
    current_gene_w_gene_set <-
        gene_w_gene_set %>%
        filter(sterility_gene_set == gene_set) 

    table <-
        genes_w_type_w_cluster %>%
        left_join(current_gene_w_gene_set) %>%
        mutate(sterility_gene=
            if_else(
                sterility_gene_set %>% is.na(),
                "Remaining Genes",
                "Sterility Genes"
            ),
            sterility_gene_set=NULL
        )

    plot <-
        ggplot(table, aes(x=cluster %>% factor(),y=sterility_gene,colour=type)) +
        geom_jitter() +
        geom_vline(xintercept = c(1:(to_keep-1)) + .5, size=.1) +
        labs(title = gene_set,
            x = "Clusters",
            y = NULL,
            colour = "Gene Type") +
        scale_colour_manual(values=pallete[c(3,1)])+
        theme_linedraw() + 
        theme(
            panel.grid.major = element_blank() 
        )

    print(plot)
}

dev.off()

