library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

path = "./"

# Input gene list
gene_list <- read.table(paste0(path, "/Sirt7.sup.table3.txt"))[,1]

# Transform to ENTREZID
gene_list.bitr <- bitr(gene_list, fromType = "SYMBOL",
                          toType = c("ENTREZID"),
                          OrgDb = org.Hs.eg.db)

birt <- gene_list.bitr

# GO enrichment
birt.ego.all <- enrichGO(
  gene = birt$ENTREZID,
  keyType = "ENTREZID",
  OrgDb = "org.Hs.eg.db",
  ont = "ALL",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05,
  readable = TRUE,
  minGSSize = 10,
  maxGSSize = 500
)

birt.ego.all.result <- birt.ego.all@result[order(birt.ego.all@result$pvalue),]
birt.ego.all.sig <- rbind(head(subset(birt.ego.all.result, ONTOLOGY == "BP"), 10),
                                       head(subset(birt.ego.all.result, ONTOLOGY == "CC"), 10),
                                       head(subset(birt.ego.all.result, ONTOLOGY == "MF"), 10))

birt.ego.all.sig$pvalue <- -log10(birt.ego.all.sig$pvalue)
birt.ego.all.sig$Description <- factor(birt.ego.all.sig$Description, levels = as.character(rev(birt.ego.all.sig$Description)))

pdf(paste0(path, "./GO.all.pdf"), width = 7, height = 6)
ggplot(birt.ego.all.sig, aes(x = Description, y = pvalue, fill = ONTOLOGY)) +
  geom_col(width = 0.5) +
  coord_flip() +
  scale_fill_manual(values = c("Red","Orange","Blue")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme_classic() +
  theme(axis.text=element_text(size=10,face="bold")) +
  facet_grid(ONTOLOGY~., scales = "free") 
dev.off()

birt.ekegg.all <- enrichKEGG(
  gene = birt$ENTREZID,
  organism  = 'hsa',
  pvalueCutoff  = 0.1,
  pAdjustMethod  = "BH",
  qvalueCutoff  = 0.1,
  minGSSize = 10,
  maxGSSize = 500
)

birt.ekegg.all.result <- as.data.frame(birt.ekegg.all)
dat <- birt.ekegg.all.result[birt.ekegg.all.result$pvalue < 0.01,]

dat$pvalue <- -log10(dat$pvalue)
dat <- dat[order(dat$pvalue, decreasing = F),]
dat <- dat[1:20,]

pdf(paste0(path, "/Users/Karl/Desktop/KEGG.pdf"),width = 7)
ggplot(dat, aes(x = reorder(Description, order(pvalue, decreasing=F)), y = pvalue, fill = pvalue)) + 
  geom_bar(stat = "identity", width = 0.5) + 
  scale_x_discrete(name = NULL) +
  scale_y_continuous(name = "-log10Pvalue") +
  scale_fill_gradientn(colours=c("pink","brown3","brown4")) +
  coord_flip() + 
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA)) +
  theme_classic() + 
  theme(plot.title = element_text(hjust = 0.5),
       axis.text=element_text(size=10,face="bold"))
dev.off()
