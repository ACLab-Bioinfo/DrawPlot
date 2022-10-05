Genelist <- c(
  paste0("Sirt",1:7)
)

seu.mouse <- readRDS("./")

for(gene in Genelist){
  message(gene)
  x <- as.vector(seu.mouse@assays$RNA@data[gene,])
  Res <- as.vector(seu.mouse$Res)
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
  x <- seu.mouse@assays$RNA@data[gene,seu.mouse@meta.data$Res == "Responder"]
  y <- seu.mouse@assays$RNA@data[gene,seu.mouse@meta.data$Res == "Non-responder"]
  mean(x)
  mean(y)
  wilcox.test(x,y)
}
