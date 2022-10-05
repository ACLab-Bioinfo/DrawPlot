############ Density UMAPPlot using ggplot
library(ggpointdensity)
library(viridis)

caf[["umap"]]
caf$anno5 <- caf$RNA_snn_res.1.1

# CAF responder
data <- cbind(Embeddings(object=caf_R[['umap']]), FetchData(caf_R, 'anno5'))
# CAF non-responder
data <- cbind(Embeddings(object=caf_NR[['umap']]), FetchData(caf_NR, 'anno5'))

p <- ggplot(data = data, mapping = aes(x = UMAP_1, y = UMAP_2)) + 
            geom_pointdensity() +
            scale_color_viridis() + theme_bw() +  
            ggtitle('NR_On_Core') + theme(plot.title = element_text(hjust = 0.5), axis.text.x = element_text(angle = 90))

p
