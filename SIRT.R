Genelist <- c(
  paste0("Sirt",1:7)
)

seu.mouse <- readRDS("/public/processed/seu.mouse-1025.rds")
seu.mouse$Res <- ifelse(seu.mouse$mouse %in% c("Mouse1", "Mouse2"), "Sensitive", "Resistant")

saveRDS(seu.mouse, "./patrick/mouse/seu.mouse-1025.051022.rds")


seu.mouse.tumor <- subset(seu.mouse, subset = anno1 == "Tumor")

# The data has not been normalized before.
seu.mouse.tumor <- NormalizeData(seu.mouse.tumor)

VlnPlot(seu.mouse.tumor, features = c("Sirt1"), split.by = "Res")

for(gene in Genelist){
  message(gene)
  x <- as.vector(seu.mouse.tumor@assays$RNA@data[gene,])
  Res <- as.vector(seu.mouse.tumor$Res)
  df <- data.frame(Res = Res, x = x)
  p <- ggplot(df,aes(x=Res, y=x, color=Res)) +
    geom_jitter(width = 0.1, size = 1, alpha = 0.8) +
    scale_color_manual(values=c("#EE3B3B", "#0000EE")) +
    theme_bw() +
    xlab("") + ylab("Log-Normalized Expression Level") +
    ggtitle(gene) +
    ylim(0,3.5) +
    theme(
      axis.text.x = element_text(size = 10, face = "plain", color = "black"),
      #    axis.text.y = element_text(size = 12, face = "plain"),
      #    axis.title = element_text(size = 12, face = "plain"),
      plot.title = element_text(hjust = 0.5))
  png(paste0("./",gene,".expression.png"), width = 5,height = 4, units='in', res=600)
  print(p)
  dev.off()
}

# statistic
for(gene in Genelist){
  x <- seu.mouse.tumor@assays$RNA@data[gene,seu.mouse.tumor@meta.data$Res == "Sensitive"]
  y <- seu.mouse.tumor@assays$RNA@data[gene,seu.mouse.tumor@meta.data$Res == "Resistant"]
  message(mean(x))
  message(mean(y))
  wil <- wilcox.test(x,y)
  print(wil)
}
