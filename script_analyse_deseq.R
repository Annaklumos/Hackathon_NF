library(DESeq2)
library(factoextra)
library(FactoMineR)

# Chargement des donnees
files <- (Sys.glob("*.counts"))
list_files <- lapply(files, function(x) read.table(x, header = T))

######  CREATION DATAFRAME
nb_ech = length(list_files)
nb_gene = nrow(list_files[[1]])

id_gene <- list_files[[1]]["Geneid"]

df <- rep(0, nb_ech)

for (i in 1:nb_ech) {
   df[i] <- data.frame(list_files[[i]][,7])
}


df <- cbind(id_gene, df)
colnames(df) <- c(1, 2, 3, 4, 5, 6, 7, 8, 9)
df_t <- t(as.matrix(df[,-1]))


#######  PCA DONNEES BRUTES #########
pdf("PCA_not_normalized.pdf")
fviz_pca_ind(PCA(df_t,graph=F),col.ind = Type)
dev.off()

#### PCA AVEC DATA FRAME AVEC UNIQUEMENT LES GENES EXPRIMES
ind_exp <- rowSums(df[,-1])
ind_exp <- which(ind_exp > 0)
df_exp <- df[ind_exp,]
id_gene_exp <- id_gene[ind_exp,]

#Analyse différentielle
Des <- DESeqDataSetFromMatrix(countData = df_exp[-1], colData = Type, design = ~Type)
analysis = DESeq(Des)
res = results(analysis)
write.table(res, "Deseq2_results_4mutants_table.txt", sep = "\t", row.names = TRUE, col.names = TRUE)

#Normalisation
vst = getVarianceStabilizedData(analysis)
vst_t <- t(vst)

#PCA SUR DONNÉES NORMALISÉES
pdf("PCA_expressed_normalized.pdf")
fviz_pca_ind(PCA(vst_t,graph=F),col.ind = Type)
dev.off()

#On voit que le mutant R625H se comporte comme un wild type, nous avons donc changé son type pour wild type
Type2 = t(data.frame("Mutant","Mutant","WT","WT","WT","WT","WT","Mutant"))
row.names(Type2) <- c(2, 3, 4, 5, 6, 7, 8, 9)
colnames(Type2) <- "Type"
Des2 <- DESeqDataSetFromMatrix(countData = df_exp[-1], colData = Type2, design = ~Type)
analysis2 = DESeq(Des2)
res2 = results(analysis2)
write.table(res2, "Deseq2_results_3mutants_table.txt", sep = "\t", row.names = TRUE, col.names = TRUE)


res_mat <- as.matrix(res2)
res_mat <- data.frame(id_gene_exp, res_mat)

### EXTRACTION DES GENES

# Lignes dont la pvalue ajustee est < 0.05
ind_padj <- res_mat[,"padj"]
ind_padj <- which(ind_padj<0.05)

exprime_diff = length(which(res_mat[,"padj"]<0.05 & (res_mat[,"log2FoldChange"]>1 |res_mat[,"log2FoldChange"]< -1)))
sur_exprime =length(which(res_mat[,"padj"]<0.05 & res_mat[,"log2FoldChange"]>1))
sous_exprime=length(which(res_mat[,"padj"]<0.05 & res_mat[,"log2FoldChange"]< -1))
print(paste("Genes sous et sur exprimés:",exprime_diff))
print(paste("Genes sous exprimés:",sous_exprime))
print(paste("Genes sur exprimés:",sur_exprime))




# Volcano plot
logp <- -log10(res_mat[,"padj"])
pdf("VolcanoPlot.pdf")
plot(logp~res_mat[,"log2FoldChange"], main = "Expression différentielle des gènes entre WT et mutés", xlab = "Log2FoldChange", ylab = "- Log10(padj)")
dev.off()