
```{r}

counts_all <- read.delim("/rds/projects/c/croftap-mapjagb9/bulkseq_JIAfibs/count/counts_all_ctrl_TGF.txt", row.names=1, comment.char="#")

library(DESeq2)

counts_all <- counts_all %>% dplyr::select(-c(Chr, Start,End, Strand, Length))
colnames(counts_all) <- c("CTRL29", "CTRL31", "CTRL25", "TGF25", "TGF29", "TGF31")



meta_data=colnames(counts_all)
meta_data<-as.data.frame(counts_all)
meta_data$conditon<-c("CTRL", "CTRL", "CTRL", "TGF", "TGF", "TGF")

dds <- DESeqDataSetFromMatrix(countData = counts_all,
                                  colData = meta_data,
                                  design = ~1)
    

gene_sub <- 200
gene_sub_max <- quantile(rowSums(counts(dds)), 0.999)

keep <- rowSums(counts(dds)) > gene_sub & rowSums(counts(dds)) < gene_sub_max
dds <- dds[keep, ]

dds@colData[['conditon']] <- as.factor(dds@colData[['conditon']])
design(dds) <- formula(~ conditon)
dds <- DESeq(dds, test = "Wald")


cond <- "conditon"

targetvar <- "TYPE"

comparisons <- data.frame(t(combn(unique(as.character(meta_data[[targetvar]])), 2)))

gene_comp <- apply(comparisons, 1, function(cp) {
  comparison <- paste0(cp[1], "_vs_", cp[2])
  degs <- as.data.frame(results(dds, contrast = c(targetvar, cp[1], cp[2])))
  degs[["gene"]] <- rownames(res)
  degs[["comparison"]] <- comparison
  degs
})

all_degs <- Reduce(rbind, gene_comp)

all_degs %>% filter(padj < 0.05) %>%
  mutate('score' = log2FoldChange*(-log10(pvalue))) %>%
  arrange(desc(abs(score))) -> all_degs_f


deseq2VST <- vst(dds, blind=T)
gene_filter <- unique(all_degs_f$gene)

deseq2VST <- assay(deseq2VST)
deseq2VST<-as.matrix(deseq2VST)
sub_vst_mat <- deseq2VST[rownames(deseq2VST) %in% gene_filter, ]
scale_sub_vst_mat <- t(scale(t(sub_vst_mat)))

colours <- list('condition' = meta_data$conditon)
col_ann <- HeatmapAnnotation(df = meta_data[,c("conditon")], col=colours)


library(colorRamp2)
Heatmap(scale_sub_vst_mat, 
              top_annotation = col_ann, 
              cluster_columns = F,
              cluster_rows = T,
              show_row_names = F,
              show_column_names = T)

table(fibs_markers_named$cluster)
fibs_markers_named_SOX5 <- fibs_markers_named %>% filter(cluster == "SOX5_CDH11")

annotation_gs <- fetchAnnotation(species="hs", ensembl_version=NULL, ensembl_host=NULL)

index <- match(subres$gene, annotation_gs$ensembl_id)
subres$gene_name <- annotation_gs$gene_name[index]

subres_f <- subres %>% filter(comparison %in% c("CTRL_vs_TGF", "CTRL_vs_EGF_TGF", "CTRL_vs_EGF") & log2FoldChange < -0.0001)


suubres_overlap <- subres_f[subres_f$gene_name %in% fibs_markers_named_SOX5$gene,]

df_overlap <- table(suubres_overlap$comparison) %>% as.data.frame()


df2 <- table(subres_f$comparison) %>% as.data.frame()
#df2 <- df2[-4,]
df_overlap$Freq.y <- df2$Freq

df_overlap$pct <- df_overlap$Freq/df_overlap$Freq.y*100


ggplot(df_overlap, aes(x=Freq, y=pct)) + 
  geom_point(size=2)+theme_ArchR()+ geom_text_repel(aes(label = df_overlap$Var1),
                    size = 3.5) +
  geom_hline(yintercept=max(df_overlap$pct)/2, linetype='dotted', col = 'red', size=1.5)+geom_vline(xintercept = max(df_overlap$Freq)/2, linetype="dotted", 
                color = "red", size=1.5)+ggtitle("DEGs with SOX5/CDH11 markers")


#pull each markers and loop through them

markers_list <- list()
for (i in 1:length(unique(fibs_markers_named$cluster))){
  
  markers_list[[i]] <- fibs_markers_named %>% filter(cluster == unique(fibs_markers_named$cluster)[[i]])

}

names(markers_list) <- unique(fibs_markers_named$cluster)

ggplots <- list()
for (i in 1:length(markers_list)){
suubres_overlap_tmp <- subres_f[subres_f$gene_name %in% markers_list[[i]]$gene,]
df_overlap_tmp <- table(suubres_overlap_tmp$comparison) %>% as.data.frame()
df2_tmp <- table(subres_f$comparison) %>% as.data.frame()
df_overlap_tmp <- merge(x = df_overlap_tmp, y = df2_tmp, by = "Var1", all = TRUE)
df_overlap_tmp <- df_overlap_tmp %>% replace(is.na(.), 0)
df_overlap_tmp$pct <- df_overlap_tmp$Freq.x/df_overlap_tmp$Freq.y*100
ggplots[[i]] <- ggplot(df_overlap_tmp, aes(x=Freq.x
, y=pct)) + 
  geom_point(size=2)+theme_ArchR()+ geom_text_repel(aes(label = df_overlap_tmp$Var1),
                    size = 3.5) +
  geom_hline(yintercept=max(df_overlap_tmp$pct)/2, linetype='dotted', col = 'red', size=1.5)+geom_vline(xintercept = max(df_overlap_tmp$Freq.x
)/2, linetype="dotted", 
                color = "red", size=1.5)+ggtitle(paste0("overlap with ", names(markers_list)[[i]]))+xlim(0,80)+ylim(0,6)



print(ggplot(df_overlap_tmp, aes(x=Freq.x
, y=pct)) + 
  geom_point(size=2)+theme_ArchR()+ geom_text_repel(aes(label = df_overlap_tmp$Var1),
                    size = 3.5) +
  geom_hline(yintercept=max(df_overlap_tmp$pct)/2, linetype='dotted', col = 'red', size=1.5)+geom_vline(xintercept = max(df_overlap_tmp$Freq.x
)/2, linetype="dotted", 
                color = "red", size=1.5)+ggtitle(paste0("overlap with ", names(markers_list)[[i]]))+xlim(0,80)+ylim(0,6))



}



plot_grid(ggplots[[2]], ggplots[[3]], ggplots[[4]], ggplots[[5]], ggplots[[6]], ggplots[[7]], ggplots[[8]], ggplots[[10]], ggplots[[11]], ncol=3, nrow=3)

subres_f$log2FoldChange %>% max

markers_overlap <- fibs_markers_named[fibs_markers_named$gene %in% subres_f$gene_name,]
markers_overlap_f <- markers_overlap %>% filter(p_val_adj < 0.05)
fibs_markers_named_f <- fibs_markers_named %>% filter(p_val_adj<0.05)
df <- table(markers_overlap_f$cluster) %>% as.data.frame()
df$total <- table(fibs_markers_named$cluster)


df$pct <- df$Freq/df$total*100

#my_patterns=c("SL_", "LL_", "LS_")
#df <-df %>%  filter(grepl(paste(my_patterns, collapse='|'), Var1))

ggplot(df, aes(x=Freq, y=pct)) + 
  geom_point()+theme_ArchR(baseSize = 12)+ geom_text_repel(aes(label = df$Var1),
                    size = 4) +
  geom_hline(yintercept=max(df$pct)/2, linetype='dotted', col = 'red', size=1)+geom_vline(xintercept = max(df$Freq
)/2, linetype="dotted", 
                color = "red", size=1)




for (i in 1:length(unique(subres_f$comparison))){
subres_f_tmp <- subres_f %>% filter(comparison == unique(subres_f$comparison)[[i]] & log2FoldChange < -2)
stroma_clean_h <- AddModuleScore(stroma_clean_h, features=list(subres_f_tmp$gene_name), name=unique(subres_f$comparison)[[i]]  )
}

stroma_clean_h$clusters
Idents(stroma_clean_h) <- 'clusters'

DotPlot(stroma_clean_h, features = c("CTRL_vs_TGF1"), idents = levels(stroma_clean_h)[c(2,3,4,6,7,8,11)])+RotatedAxis()

dotplot<-DotPlot(stroma_clean_h, features = c("CTRL_vs_TGF1", "TGFBR1", "TGFBR2", "TGFBR3"), idents = levels(stroma_clean_h)[c(2,3,4,6,7,8,11)])+RotatedAxis()

dotplot<-dotplot$data

dotplot<-dotplot %>% 
  dplyr::select(-pct.exp, -avg.exp) %>%  
  pivot_wider(names_from = id, values_from = avg.exp.scaled) %>% 
  as.data.frame() 

dotplot$features.plot<-unique(dotplot$features.plot)
dotplot<-na.omit(dotplot)

row.names(dotplot) <- dotplot$features.plot  
dotplot <- dotplot[,-1] %>% as.matrix()

library(ComplexHeatmap)

Heatmap(dotplot, border=T)




























