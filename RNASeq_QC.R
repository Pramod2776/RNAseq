
require(sva)
require(ggplot2)
require(edgeR)
require(reshape2))
library("ggfortify")
library(tidyverse)  
library(cluster)  
library(factoextra) 
library(dendextend)
library(diceR)
library(dplyr)
library(ggplot2)
library(sigclust)
library(rtracklayer)
library(data.table)
library(ComplexHeatmap)
library(circlize)

## Sample meta data  and RSEM

sampleinfo <- read.delim("./Sample_Metadata.txt", 
                         header = T, check.names = F)


counts <- read.delim("./rsem.merged.gene_counts.tsv",
                     header = T, check.names = F, row.names = 1)

gene2tx <- data.frame(`gene_id`=rownames(counts), 
                      `transcript_id(s)`=counts$`transcript_id(s)`)
counts <- counts[,2:dim(counts)[2]]


## Boxplot and PCA plot of Raw data

tpm <- read.delim("./rsem.merged.gene_tpm.tsv",
                  header = T, check.names = F, row.names = 1)
tpm <- tpm[,2:dim(tpm)[2]]


## Filter lowly expressed genes from TPM 


tpm.cutoff <- ncol(tpm) * 0.1
tpm.sel <- tpm[apply(tpm,1,sum) >= tpm.cutoff, ]
tpm.sel.log2 <- log2(tpm.sel + 1.0)
tpm.sel.log2.melt <- melt(as.matrix(tpm.sel.log2))
tpm.sel.log2.melt$Batch <- sampleinfo$Batch[match(tpm.sel.log2.melt$Var2, sampleinfo$Sample)]
head(tpm.sel.log2.melt)

p1 <- ggplot(tpm.sel.log2.melt, aes(x=Var2, y=value)) + 
  geom_boxplot(aes(fill=Batch, alpha=0.90), outlier.shape=NA)  + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_fill_manual(values=c("red2","blue2")) + xlab("Sample") + labs(title="") +
  ylab("Log2 (TPM)") +
  scale_alpha(guide = 'none') + ylim(0,7.5)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# PCA plot
expr.pca <- t(tpm.sel.log2)
expr.pca <- expr.pca[, apply(expr.pca, 2, var, na.rm = T) != 0]
expr.pca <- prcomp(expr.pca, center = T, scale. = T)
Batch <- sampleinfo$Batch
pca.groups <- data.frame(Batch)
head(pca.groups)
pca.groups <- cbind(pca.groups, as.data.frame(expr.pca$x[,1:3]))

cols <- c("red2","blue2")
autoplot(expr.pca, x=1, y=2, data=pca.groups, colour="Batch", size = 4, alpha = 0.90) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 

autoplot(expr.pca, x=1, y=3, data=pca.groups, colour="Batch", size = 4, alpha = 0.90) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 


# batch effect Correction

adjusted_counts <- ComBat_seq(as.matrix(counts), batch=Batch)

# Filter genes
count.cutoff <- ncol(adjusted_counts) * 5
adjusted_counts.sel <- adjusted_counts[apply(adjusted_counts,1,sum) >= count.cutoff, ]
dim(adjusted_counts.sel)
gene2tx.sel <- gene2tx[match(rownames(adjusted_counts.sel), gene2tx$gene_id),]
table(gene2tx.sel$gene_id == rownames(adjusted_counts.sel))
output1 <- cbind(gene2tx.sel, adjusted_counts.sel)

# calculate RPKM
rsem_output <- read.delim("../sample.genes.results",
                          header = T)

gene.lengths <- rsem_output$length[match(rownames(adjusted_counts), rsem_output$gene_id)]
head(gene.lengths)

rpkms <- rpkm(adjusted_counts, gene.length = gene.lengths)
rpkms <- round(rpkms, digits = 2)


table(rownames(rpkms) == gene2tx$gene_id)
outputa <- cbind(gene2tx, rpkms)

rpkms.sel <- rpkms[rownames(adjusted_counts.sel),]

table(gene2tx.sel$gene_id == rownames(rpkms.sel))
output4 <- cbind(gene2tx.sel, rpkms.sel)

# Boxplot 

rpkms.sel.log2 <- log2(rpkms.sel + 1.0)
rpkms.sel.log2.melt <- melt(rpkms.sel.log2)

rpkms.sel.log2.melt$Batch <- sampleinfo$Batch[match(rpkms.sel.log2.melt$Var2, sampleinfo$Sample)]

p1 <- ggplot(rpkms.sel.log2.melt, aes(x=Var2, y=value)) + 
  geom_boxplot(aes(fill=Batch, alpha=0.90), outlier.shape=NA)  + theme_classic() + 
  theme(panel.border = element_rect(colour = "black", fill=NA, size=1))+
  scale_fill_manual(values=c("red2","blue2")) + xlab("Sample") + labs(title="") +
  ylab("Log2 (RPKM + 1.0)") +
  scale_alpha(guide = 'none') + ylim(0,10)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
print(p1)


# PCA plot
expr.pca <- t(rpkms.sel.log2)
expr.pca <- expr.pca[, apply(expr.pca, 2, var, na.rm = T) != 0]
expr.pca <- prcomp(expr.pca, center = T, scale. = T)

pca.groups <- data.frame(Batch)
head(pca.groups)
pca.groups <- cbind(pca.groups, as.data.frame(expr.pca$x[,1:3]))

cols <- c("red2","blue2")

autoplot(expr.pca, x=1, y=2, data=pca.groups, colour="Batch", size = 4, alpha = 0.90) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 

autoplot(expr.pca, x=1, y=2, data=pca.groups, colour="Batch", shape = FALSE, size = 4, alpha = 0.90, label.size = 5) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 

autoplot(expr.pca, x=1, y=3, data=pca.groups, colour="Batch", size = 4, alpha = 0.90) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 

autoplot(expr.pca, x=1, y=3, data=pca.groups, colour="Batch", shape = FALSE, size = 4, alpha = 0.90, label.size = 6) + 
  theme_classic() + theme(panel.border = element_rect(colour = "black", fill=NA, size=1)) +
  scale_color_manual(values = cols) 


# clustering of samples using batch corrected RPKM #

df <- scale(t(rpkms.sel.log2))
d <- dist(df, method = "euclidean")
hc1 <- hclust(d, method = "complete" )

# Plot the obtained dendrogram
plot(hc1, cex = 0.6, hang = -1)

# Compute with agnes
hc2 <- agnes(df, method = "complete")

# Agglomerative coefficient
hc2$ac
## [1] 0.5770707

# methods to assess
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

# function to compute coefficient
ac <- function(x) {
  agnes(df, method = x)$ac
}

map_dbl(m, ac)
##average    single  complete      ward 
##  0.5257359 0.4921319 0.5770707 0.7962639 

hc3 <- agnes(df, method = "ward")
pltree(hc3, cex = 0.6, hang = -1, main = "Dendrogram of agnes")

# compute divisive hierarchical clustering
hc4 <- diana(df)

# Divise coefficient; amount of clustering structure found
hc4$dc
## [1] 0.5631029

# plot dendrogram
pltree(hc4, cex = 0.6, hang = -1, main = "Dendrogram of diana")

# Ward's method
hc5 <- hclust(d, method = "ward.D2" )

# Cut tree into 4 groups
sub_grp <- cutree(hc5, k = 10)

# Number of members in each cluster
table(sub_grp)
## sub_grp
##  1  2  3  4 
##  17 20 8  3

clust_df = df %>%
  as.data.frame() %>%
  mutate(cluster = sub_grp) %>%
  select(cluster, everything()) 

plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 10, border = 2:5)

pdf(file = "GBM_hclust_ward.D2.pdf", height = 6, width = 7)
plot(hc5, cex = 0.6)
rect.hclust(hc5, k = 10, border = 2:5)
dev.off()

pdf(file = "GBM_hclust_group_ward.D2.pdf", height = 8, width = 12)
fviz_cluster(list(data = df, cluster = sub_grp))
dev.off()


# Cut agnes() tree into 4 groups
hc_a <- agnes(df, method = "ward")
cutree(as.hclust(hc_a), k = 4)

# Cut diana() tree into 4 groups
hc_d <- diana(df)
cutree(as.hclust(hc_d), k = 4)

# Compute distance matrix
res.dist <- dist(df, method = "euclidean")

# Compute 2 hierarchical clusterings
hc1 <- hclust(res.dist, method = "complete")
hc2 <- hclust(res.dist, method = "ward.D2")

# Create two dendrograms
dend1 <- as.dendrogram (hc1)
dend2 <- as.dendrogram (hc2)

tanglegram(dend1, dend2)

dend_list <- dendlist(dend1, dend2)

tanglegram(dend1, dend2,
           highlight_distinct_edges = FALSE, # Turn-off dashed lines
           common_subtrees_color_lines = FALSE, # Turn-off line colors
           common_subtrees_color_branches = TRUE, # Color common branches 
           main = paste("entanglement =", round(entanglement(dend_list), 2))
)


fviz_nbclust(df, FUN = hcut, method = "wss")

fviz_nbclust(df, FUN = hcut, method = "silhouette")

gap_stat <- clusGap(df, FUN = hcut, nstart = 25, K.max = 10, B = 50)
fviz_gap_stat(gap_stat)

## heatmap

mat = na.omit(mat)
mat$Gene_name = NULL
mat$Ensembl_ID = NULL
mat$Gene_biotype = NULL
colnames(mat)

mat = mat %>%
  as.data.frame() %>%
  dplyr::select(Gene_name, Subtype, everything())


mat_score =  mat[,3:50]
rownames(mat_score) = mat$Gene_name
base_mean = rowMeans(mat_score)
mat_scaled = t(apply(mat_score, 1, scale))

max(mat_scaled)
min(mat_scaled)

mat_scaled [mat_scaled > 2.0] <- 2.0
mat_scaled [mat_scaled < -2.0] <- -2.0 

colors <- colorRamp2(c(-2, 0, 2), c("royalblue","white", "red3"))

table(mat$Subtype)

Category <- factor(mat$Subtype, levels = c("C", "M", "N", "P", "U"))

colnames(mat[,3:50])
colnames(mat_scaled) = colnames(mat[,3:50])

pdf(file = "Heatmap_mat.pdf", height = 10,width = 10)
ht1 = Heatmap(mat_scaled, col = colors, row_split = Category,
              border = T,na_col = "#737373",cluster_columns = T,
              show_row_names = F, show_column_names =T,cluster_rows = T,
              row_dend_side = "left",column_names_rot = 90,
              column_dend_side = "top",
              row_names_max_width = unit(9, "cm"),
              row_gap = unit(2, "mm"),
              heatmap_legend_param = list(at = c(-2,-1,0,1,2)),
              row_dend_width = unit(1.5, "cm")) 

ht1
dev.off()





