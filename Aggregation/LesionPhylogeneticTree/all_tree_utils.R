library(ape)      # for phylo trees
library(igraph)  # for graphs

pretty_tree <- function(dataset, clusters, dist_method = "euclidean", clus_method = "complete") {
    # distance matrix
    ## hcluster = hclust(dist_data, method = clus_method)
    dist_data = dist(dataset, method = dist_method)
    # hierarchical clustering
    hcluster = hclust(dist_data, method = clus_method)
    
    # convert to phylo object
    phylo_tree = as.phylo(hcluster)
    # get edges
    graph_edges = phylo_tree$edge
    # convert to graph
    graph_net = graph.edgelist(graph_edges)
    # extract layout (x-y coords)
    graph_layout = layout.auto(graph_net)
    
    # colors like default ggplot2
    ggcolors <- function(n, alfa) {
        hues = seq(15, 375, length = n + 1)
        hcl(h = hues, l = 65, c = 100, alpha = alfa)[1:n]
    }
    
    # colors of labels and points
    num_clusters <- length(unique(clusters))
    txt_pal = ggcolors(num_clusters)
    pch_pal = paste(txt_pal, "55", sep='')
    txt_col = txt_pal[clusters]
    pch_col = pch_pal[clusters]
    
    # additional params
    nobs = length(clusters)
    nedges = nrow(graph_edges)
    
    # start plot
    plot(graph_layout[,1], graph_layout[,2], type = "n", axes = FALSE,
         xlab = "", ylab = "")
    # draw tree branches
    segments(
        x0 = graph_layout[graph_edges[,1],1],
        y0 = graph_layout[graph_edges[,1],2],
        x1 = graph_layout[graph_edges[,2],1],
        y1 = graph_layout[graph_edges[,2],2],
        col = "black", lwd = 5
    )
    # add tree leafs
    points(graph_layout[1:nobs,1], graph_layout[1:nobs,2], col = pch_col,
           pch = 19, cex = 2)
    # add empty nodes
    points(graph_layout[(nobs+1):nedges,1], graph_layout[(nobs+1):nedges,2],
           col = "black", pch = 19, cex = 0.5)
}