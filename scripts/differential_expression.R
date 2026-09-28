suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
})
set.seed(73681)
dir.create("results/differential_expression", recursive=TRUE, showWarnings=FALSE)
dir.create("results/figures", recursive=TRUE, showWarnings=FALSE)
meta <- read.delim("metadata/samples.tsv", check.names=FALSE, stringsAsFactors=FALSE)
cts <- as.matrix(read.delim("results/counts/count_matrix.tsv", row.names=1, check.names=FALSE))
stopifnot(identical(colnames(cts), meta$sample), !anyDuplicated(rownames(cts)),
          all(is.finite(cts)), all(cts>=0), all(cts==round(cts)), ncol(cts)==6)
storage.mode(cts) <- "integer"
rownames(meta) <- meta$sample
meta$condition <- factor(meta$condition, levels=c("WT", "sec66del"))
stopifnot(all(table(meta$condition)==3))
# Replicate numbers are independent isolates, not evidence of pairing or a batch variable.
dds <- DESeqDataSetFromMatrix(countData=cts, colData=meta, design=~condition)
keep <- rowSums(counts(dds)>=10)>=3
dds <- DESeq(dds[keep,])
res <- results(dds, contrast=c("condition", "sec66del", "WT"), alpha=0.05)
out <- as.data.frame(res)
out$gene_id <- rownames(out)
ann <- read.delim("metadata/gene_annotation.tsv", stringsAsFactors=FALSE)
out$gene_name <- ann$gene_name[match(out$gene_id, ann$gene_id)]
out$biotype <- ann$biotype[match(out$gene_id, ann$gene_id)]
out$significant_FDR05 <- !is.na(out$padj) & out$padj<0.05
# Additional effect-size display filter; the Wald test itself tests zero log2 fold change.
out$large_effect_FDR05 <- out$significant_FDR05 & abs(out$log2FoldChange)>=1
out <- out[order(out$padj, na.last=TRUE),c("gene_id","gene_name","biotype",names(out)[1:6],
                                       "significant_FDR05","large_effect_FDR05")]
write.table(out, "results/differential_expression/all_tested_genes.tsv", sep="\t", quote=FALSE,row.names=FALSE)
write.table(out[out$large_effect_FDR05,], "results/differential_expression/DEGs_FDR05_absLFC1.tsv",
            sep="\t", quote=FALSE,row.names=FALSE)
write.table(data.frame(gene_id=rownames(dds),counts(dds,normalized=TRUE),check.names=FALSE),
            "results/counts/normalized_counts.tsv", sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(sample=colnames(dds),size_factor=sizeFactors(dds)),
            "results/differential_expression/size_factors.tsv", sep="\t",quote=FALSE,row.names=FALSE)

vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
pca <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percent <- round(100*attr(pca,"percentVar"),1)
write.table(pca,"results/differential_expression/PCA_coordinates.tsv",sep="\t",quote=FALSE,row.names=FALSE)
p <- ggplot(pca,aes(PC1,PC2,color=condition,label=name)) + geom_point(size=3.5) +
  geom_text(vjust=-0.8,size=3.2,show.legend=FALSE) +
  scale_color_manual(values=c(WT="#2878A0",sec66del="#D26742")) +
  scale_x_continuous(expand=expansion(mult=c(0.12,0.18))) +
  scale_y_continuous(expand=expansion(mult=c(0.10,0.16))) +
  labs(title="SEC66 deletion: sample relationships",subtitle="GSE73681 · six complete RNA-seq libraries",
       x=paste0("PC1 (",percent[1],"%)"), y=paste0("PC2 (",percent[2],"%)"),color="Condition") +
  theme_minimal(base_size=12) + theme(legend.position="bottom",plot.margin=margin(15,25,15,15))
ggsave("results/figures/PCA.png",p,width=8,height=5.5,dpi=180)

sample_distances <- as.matrix(dist(t(assay(vsd))))
write.table(sample_distances,"results/differential_expression/sample_distances.tsv",
            sep="\t",quote=FALSE,col.names=NA)
distance_tree <- hclust(as.dist(sample_distances))
pheatmap(sample_distances,cluster_rows=distance_tree,cluster_cols=distance_tree,
         main="Sample distances · variance-stabilised counts",
         filename="results/figures/sample_distances.png",width=7,height=6)

png("results/figures/MA.png",width=1440,height=990,res=180)
plotMA(res,alpha=0.05,main="sec66 deletion vs wild type (FDR < 0.05)",ylim=c(-6,6))
dev.off()
v <- out[!is.na(out$padj) & is.finite(out$log2FoldChange),]
v$status <- "Other tested genes"
v$status[v$large_effect_FDR05 & v$log2FoldChange>0] <- "Higher in sec66 deletion"
v$status[v$large_effect_FDR05 & v$log2FoldChange<0] <- "Lower in sec66 deletion"
v$neglog10padj <- -log10(pmax(v$padj,1e-300))
p <- ggplot(v,aes(log2FoldChange,neglog10padj,color=status))+geom_point(alpha=0.65,size=1.1)+
  geom_vline(xintercept=c(-1,1),linetype=2,color="grey65")+
  geom_hline(yintercept=-log10(0.05),linetype=2,color="grey65")+
  scale_color_manual(values=c("Higher in sec66 deletion"="#D26742","Lower in sec66 deletion"="#2878A0",
                              "Other tested genes"="#B3B3B3"))+
  labs(title="Differential expression",subtitle="DESeq2 · sec66 deletion relative to WT",
       x="Log2 fold change (unshrunken)",y="−log10 adjusted p-value",color=NULL)+
  theme_minimal(base_size=12)+theme(legend.position="bottom",legend.text=element_text(size=9))
ggsave("results/figures/volcano.png",p,width=8,height=5.5,dpi=180)

top <- head(out$gene_id[out$significant_FDR05],20)
if(length(top)>=2) {
  mat <- assay(vsd)[top,,drop=FALSE]
  rownames(mat) <- paste(out$gene_name[match(top,out$gene_id)],top,sep=" | ")
  annotation <- data.frame(condition=meta$condition,row.names=meta$sample)
  pheatmap(mat,scale="row",annotation_col=annotation,
           annotation_colors=list(condition=c(WT="#2878A0",sec66del="#D26742")),
           fontsize_row=8,main="Top 20 genes by adjusted p-value · row-scaled VST",
           filename="results/figures/top20_heatmap.png",width=8,height=8)
}
stats <- data.frame(metric=c("genes_in_matrix","genes_after_low_count_filter","genes_with_nonmissing_padj",
                              "FDR05_genes","FDR05_absLFC1_genes","higher_in_sec66del_absLFC1",
                              "lower_in_sec66del_absLFC1","PC1_percent","PC2_percent"),
                    value=c(nrow(cts),sum(keep),sum(!is.na(out$padj)),sum(out$significant_FDR05),
                            sum(out$large_effect_FDR05),
                            sum(out$large_effect_FDR05 & out$log2FoldChange>0,na.rm=TRUE),
                            sum(out$large_effect_FDR05 & out$log2FoldChange<0,na.rm=TRUE),percent))
write.table(stats,"results/analysis_summary.tsv",sep="\t",quote=FALSE,row.names=FALSE)
write.table(out[out$gene_id=="YBR171W",],"results/differential_expression/SEC66_check.tsv",sep="\t",quote=FALSE,row.names=FALSE)
saveRDS(dds,"results/differential_expression/dds.rds")
capture.output(sessionInfo(),file="metadata/R_sessionInfo.txt")
print(stats,row.names=FALSE)
