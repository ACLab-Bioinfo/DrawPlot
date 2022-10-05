############ Density UMAPPlot using ggplot
library(ggpointdensity)
library(viridis)

caf[["umap"]]
caf$anno5 <- caf$RNA_snn_res.1.1

# CAF responder
data <- cbind(Embeddings(object=caf_R[['umap']]), FetchData(caf_R, 'anno5'))
# CAF non-responder
data <- cbind(Embeddings(object=caf_NR[['umap']]), FetchData(caf_NR, 'anno5'))

# CAF Pre-treatment
data <- cbind(Embeddings(object=caf_Pre[['umap']]), FetchData(caf_Pre, 'anno5'))
# CAF On-treatment
data <- cbind(Embeddings(object=caf_On[['umap']]), FetchData(caf_On, 'anno5'))

# CAF Edge
data <- cbind(Embeddings(object=caf_Edge[['umap']]), FetchData(caf_Edge, 'anno5'))
# CAF Core
data <- cbind(Embeddings(object=caf_Core[['umap']]), FetchData(caf_Core, 'anno5'))

# CAF responder & Pre/On
data <- cbind(Embeddings(object=caf_NR_Pre_Core[['umap']]), FetchData(caf_NR_Pre_Core, 'anno5'))
data <- cbind(Embeddings(object=caf_NR_On_Core[['umap']]), FetchData(caf_NR_On_Core, 'anno5'))


p <- ggplot(data = data, mapping = aes(x = UMAP_1, y = UMAP_2)) + 
            geom_pointdensity() +
            scale_color_viridis() + theme_bw() +  
            ggtitle('NR_On_Core') + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90))

p
