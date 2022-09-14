library(ggplot2)
library(dplyr)

DrawTIDE <- function(file, gene){
### TIDE score
dir.create(file.path("./output", gene))
pdf(paste0("./output/", gene, "/Vlnplot.TIDEscore.pdf"), width = 4, height = 4)
p1 <- ggplot(file, aes(x=expr_high_or_low, y=TIDE, fill=expr_high_or_low)) + # y = TIDE/ Dysfunction/ Exclusion/ MDSC score
  geom_violin(width=0.7, trim = FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(name = "Group", values=c(
    c("#EE3B3B", "#0000EE")
  )) +
  theme_classic() +
  xlab("") + ylab("TIDE Score") +
  theme(
    axis.text.x = element_text(size = 15, face = "bold"), # high, low
    axis.text.y = element_text(size = 16, face = "bold"), # 2, 1, 0, -1, -2
    axis.title = element_text(size = 15, face = "bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "dark grey")
print(p1)
dev.off()

### Dysfunction score
pdf(paste0("./output/", gene, "/Vlnplot.Dysfunctionscore.pdf"), width = 4, height = 4)
p2 <- ggplot(file, aes(x=expr_high_or_low, y=Dysfunction, fill=expr_high_or_low)) + # y = TIDE/ Dysfunction/ Exclusion/ MDSC score
  geom_violin(width=0.7, trim = FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(name = "Group", values=c(
    c("#EE3B3B", "#0000EE")
  )) +
  theme_classic() +
  xlab("") + ylab("Dysfunction Score") +
  theme(
    axis.text.x = element_text(size = 15, face = "bold"), # high, low
    axis.text.y = element_text(size = 16, face = "bold"), # 2, 1, 0, -1, -2
    axis.title = element_text(size = 15, face = "bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "dark grey")
print(p2)
dev.off()

### Exclusion score
pdf(paste0("./output/", gene, "/Vlnplot.Exclusionscore.pdf"), width = 4, height = 4)
p3 <- ggplot(file, aes(x=expr_high_or_low, y=Exclusion, fill=expr_high_or_low)) + # y = TIDE/ Dysfunction/ Exclusion/ MDSC score
  geom_violin(width=0.7, trim = FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(name = "Group", values=c(c("#EE3B3B", "#0000EE"))) +
  theme_classic() +
  xlab("") + ylab("Exclusion Score") +
  theme(
    axis.text.x = element_text(size = 15, face = "bold"), # high, low
    axis.text.y = element_text(size = 16, face = "bold"), # 2, 1, 0, -1, -2
    axis.title = element_text(size = 15, face = "bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "dark grey")
print(p3)
dev.off()

### MDSC score
pdf(paste0("./output/",gene,"/Vlnplot.MDSCscore.pdf"), width = 4, height = 4)
p4 <- ggplot(file, aes(x=expr_high_or_low, y=MDSC, fill=expr_high_or_low)) + # y = TIDE/ Dysfunction/ Exclusion/ MDSC score
  geom_violin(width=0.7, trim = FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  scale_fill_manual(name = "Group", values=c("#EE3B3B", "#0000EE")) +
  theme_classic() +
  xlab("") + ylab("MDSC Score") +
  theme(
    axis.text.x = element_text(size = 15, face = "bold"), # high, low
    axis.text.y = element_text(size = 16, face = "bold"), # 2, 1, 0, -1, -2
    axis.title = element_text(size = 15, face = "bold")) +
  geom_hline(yintercept=0, linetype="dashed", color = "dark grey")
print(p4)
dev.off()

# Waterfall plot
col <- c("#EE3B3B", "#0000EE")
x <- 1:nrow(file)
pdf(paste0("./output/",gene,"/Waterfullplit.TDEscore.pdf"), width = 8, height = 4)
p5 <- ggplot(file, aes(x = x, y=TIDE, fill=expr_high_or_low, color=expr_high_or_low)) +
  geom_bar(stat="identity", width=0.4, position = position_dodge(width=0.4)) +
  scale_fill_manual(name="Group", values = col) + scale_color_manual(guide="none",values=col) +
  theme_classic() +
  xlab("") + ylab("TIDE score") +
  geom_vline(xintercept = which(file$TIDE<0)[1] - 0.5, linetype="dashed", color = "black") +
  theme(
    axis.line.x = element_blank(), axis.text.x = element_blank(),axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 16, face = "bold"), # 2, 1, 0, -1, -2
    axis.title = element_text(size = 15, face = "bold"))+
  coord_cartesian(ylim = c(-2,2)) +
  scale_x_continuous(expand = c(0, 0.5)) + 
  scale_y_continuous(expand = c(0, 0))
print(p5)
dev.off()


TIDE.x = file$TIDE[file$expr_high_or_low == "high"]
TIDE.y = file$TIDE[file$expr_high_or_low == "low"]
TIDE.stat <- t.test(TIDE.x, TIDE.y, alternative = "greater")
Dysfunction.x = file$Dysfunction[file$expr_high_or_low == "high"]
Dysfunction.y = file$Dysfunction[file$expr_high_or_low == "low"]
Dysfunction.stat <- t.test(Dysfunction.x, Dysfunction.y, alternative = "greater")
Exclusion.x = file$Exclusion[file$expr_high_or_low == "high"]
Exclusion.y = file$Exclusion[file$expr_high_or_low == "low"]
Exclusion.stat <- t.test(Exclusion.x, Exclusion.y, alternative = "greater")
MDSC.x = file$MDSC[file$expr_high_or_low == "high"]
MDSC.y = file$MDSC[file$expr_high_or_low == "low"]
MDSC.stat <- t.test(MDSC.x, MDSC.y, alternative = "greater")

H_R = file %>% 
  filter(Responder == "True" & expr_high_or_low == 'high')

H_NR = file %>% 
  filter(Responder == "False" & expr_high_or_low == 'high')

L_R = file %>% 
  filter(Responder == "True" & expr_high_or_low == 'low')

L_NR = file %>% 
  filter(Responder == "False" & expr_high_or_low == 'low')

TIDE.Response.stat <- chisq.test(data.frame(
  Responder = c(nrow(H_R),nrow(L_R)),
  Non.Responder = c(nrow(H_NR),nrow(L_NR)),
  row.names = c("High","Low"))
)

print(paste0("P-value (t test) of TIDE score is ", TIDE.stat$p.value))
print(paste0("P-value (t test) of Dysfunction score is ", Dysfunction.stat$p.value))
print(paste0("P-value (t test) of Exclusion score is ", Exclusion.stat$p.value))
print(paste0("P-value (t test) of MDSC score is ", MDSC.stat$p.value))
print(paste0("P-value (Chi-squared test) of MDSC score is ", TIDE.Response.stat$p.value))

}

file = read.csv("/Users/Karl/Downloads/R/Chengpeng/FRA1_TCGA_TIDE/first_last_25_TIDE.result.withFOSL1.expr.txt")
gene = "FRA1"
DrawTIDE(file, gene)

