########### Pearson Correlation test

# Plot Pearson Correlation between CytoSig Zscore(cytokine activity) and gene expression

# Cytokine activity .Zscore file
df <- read.table("/home/patrick/patrick/tumor/.Zscore", header = TRUE, sep = "\t")

# Gene barcode
colnames(df)


# Cytokine name
rownames(df)

# Tumor data subset from seu.0928.rds
tumor <- subset(seu, subset = anno1 == "Tumor")

# Matching of barcode between .Zscore file and tumor data
colnames(tumor) == colnames(df)

colnames(df) <- str_replace(colnames(df), "\\.", "-")

for (i in 1:43) {

# Extract one cytokine activity from the .Zscore dataframe
cytokine <- unlist(df[i,])
cytokine_name <- rownames(df[i,])

# Extract FOSL1 gene expression value from tumor data
gene <- tumor@assays$RNA@data["FOSL1",colnames(df)]

# Run the pearson correlation test
# cor.test(my_data_without_zero$gene, my_data_without_zero$cytokine, method= "pearson")

my_data <- data.frame(gene, cytokine)
my_data_without_zero <- my_data %>% 
  filter(gene!=0)

# x = gene expression (FOSL1), y = cytokine Zscore (activity)
p <- ggscatter(my_data_without_zero, x = "gene", y = "cytokine", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "FOSL1 expression level", ylab = paste0(cytokine_name," Activity"))
system(paste0("mkdir -p Cytosig_FOSL1"))
pdf(paste0('./Cytosig_FOSL1/', cytokine_name, "_cor.pdf"), width = 5, height = 5)
print(p)
dev.off()

}
