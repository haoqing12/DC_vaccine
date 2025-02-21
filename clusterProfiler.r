library(org.Hs.eg.db)
library(clusterProfiler)
# library(DOSE)
library(enrichplot)



for (cell_type in scRNA_sub_T@meta.data$Tcelltype %>% unique){

tmp <- FindMarkers(scRNA_sub_T, ident.1 = 'Post', ident.2 = 'Pre', group.by = 'status',
                       subset.ident = cell_type) %>%
      filter(p_val_adj < 0.05)

diff_entrez <- bitr(rownames(tmp),
                    fromType = "SYMBOL",
                    toType = "ENTREZID",
                    OrgDb = "org.Hs.eg.db")

go_enrich_results <- enrichGO(gene = diff_entrez$ENTREZID,
                              OrgDb = "org.Hs.eg.db",
                              ont   = "ALL"  ,     
                              readable  = TRUE)

save(go_enrich_results, file =paste0('GO/',cell_type,'_venn_GO.Rdata'))
write.csv(go_enrich_results@result, paste0('GO/',cell_type,'_venn_GO.csv'))
}