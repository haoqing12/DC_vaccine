library(data.table)
library(tidyverse)
library(maftools)


MAF <- fread("pass_maf_data.tsv")

info <- fread('/work/haoq/Neo_scRNA/clincal.tsv') 

MAF$Variant_Classification %>% table

MAF <- MAF %>% distinct %>%
        mutate(consequence = case_when(
            Variant_Classification == "Frame_Shift_Del" ~ "Frameshift",
            Variant_Classification == "Frame_Shift_Ins" ~ "Frameshift",
            Variant_Classification == "In_Frame_Del" ~ "Inframe",
            Variant_Classification == "In_Frame_Ins" ~ "Inframe",
            Variant_Classification == "Missense_Mutation" ~ "Missense",
            Variant_Classification == "Nonsense_Mutation" ~ "Nonsense",
            Variant_Classification == "Nonstop_Mutation" ~ "Nonsense",
            Variant_Classification == "Splice_Site" ~ "Splicing"
        ))

MAF$Tumor_Sample_Barcode %>% table %>% as.data.frame()

# RNA
gene_list <- fread("pass_RNA_gene.summary.tsv") %>% head(25)
summary <- fread("pass_RNA_variant.classification.summary.tsv") %>% left_join(info,by = "Tumor_Sample_Barcode")
fraction <- fread("pass_RNA_fraction.contribution.tsv")

# DNA
gene_list <- fread("pass_gene.summary.tsv") %>% head(25)
summary <- fread("pass_variant.classification.summary.tsv") %>% left_join(info,by = "Tumor_Sample_Barcode")
fraction <- fread("pass_fraction.contribution.tsv")

df0 <- data.frame()
for (y in info$Tumor_Sample_Barcode){
    d2 <- data.frame()
    for (i in gene_list$Hugo_Symbol){
        print(i)
        tmp <- MAF %>% filter(Hugo_Symbol == i & Tumor_Sample_Barcode == y)%>%
            group_by(Tumor_Sample_Barcode, Hugo_Symbol) %>%
            summarise(class = paste(consequence, collapse = ","))
            # mutate(class = ifelse(grepl(',',var), "Multi_Hit",var)) %>% ungroup()

        d1 <- data.table(gene = tmp$class) 
        # rownames(d1) <- y
        colnames(d1) <- i

    d2 <- cbind(d1,d2)
    }

    d2$sample <- y
    df0 <- rbind(df0,d2)
}


df0[is.na(df0)] = ""

df0 <- df0 %>% left_join(info, by = c("sample"="Tumor_Sample_Barcode"))

df0 %>% fwrite("pass_RNA_maf_plot.tsv", sep="\t")


#### 

# DNA
gene_list <- fread("pass_gene.summary.tsv") %>% head(25)
summary <- fread("pass_variant.classification.summary.tsv") %>% left_join(info,by = "Tumor_Sample_Barcode")
fraction <- fread("pass_fraction.contribution.tsv")

df0 <- 
  fread("/work/haoq/Neo_scRNA/pass_maf_plot.tsv") %>%  ### DNA
  mutate(phase_ = ifelse(sample %in% c("yangChunji","huangHuade","fanGuanghou","huangNian","liuDingping","yanDingbin","heQuan","tangYuanbi","wuFulong","pengJianzhong","wenJianpingII","dengShuanglongII"), "I","II"),
      phase = paste0(phase_,"_",Resp))


df1 <- df0 %>% select(TP53:RYR2, CDKN2A, PIK3CA, RYR1) %>% t()
dim(df1)
colnames(df1) <- df0$ID

col = c(
  'Missense' = "#ECB476",
  'Nonsense' = "#e14522",
  'Frameshift' = "#2cb2ca",
  'Inframe' = "#b96319",
  'Splicing' =  "#FFC001",
  'Multi_Hit' = "black"
)

alter_fun = list(
    background = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = "#CCCCCC", col = NA))
    },
    Missense = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Missense"], col = NA))
    },
     Nonsense = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Nonsense"], col = NA))
    },
    Frameshift = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Frameshift"], col = NA))
    },
    Inframe = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Inframe"], col = NA))
    },
    Splicing = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Splicing"], col = NA))
    },
    Multi_Hit = function(x, y, w, h) {
        grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
            gp = gpar(fill = col["Multi_Hit"], col = NA))
    }
    
)


col = c(
  'Missense' = "#ECB476",
  'Nonsense' = "#e14522",
  'Frameshift' = "#2cb2ca",
  'Inframe' = "#b96319",
  'Splicing' =  "#FFC001",
  'Multi_Hit' = "black"
)

library(ComplexHeatmap)

Top_column = HeatmapAnnotation(
      All_mut = anno_barplot(summary %>% select(Frame_Shift_Del:Nonstop_Mutation),bar_width = 0.7,
            gp = gpar(fill = c("#2cb2ca","#2cb2ca","#b96319","#b96319","#ECB476","#e14522","#e14522","#FFC001"), 
                    col = c("#2cb2ca","#2cb2ca","#b96319","#b96319","#ECB476","#e14522","#e14522","#FFC001"), alpha = 1),
            height = unit(1, "cm"),border = FALSE),
        # phase = df0$phase,
        Resp = df0$Resp,
        Sex = df0$Sex,
        Age = df0$Age2,
        # ECGO = df0$ECGO,
        TRS = df0$TRS,
        # ypStage = df0$ypStage,
        TG = df0$IP,
        # ID = df0$ID,
        annotation_height = c(1.5,0.2,0.2,0.2,0.2,0.2),
        # annotation_width = c(1,1,1,1,1,1),
        annotation_name_side = "left",
        annotation_name_gp = gpar(fontsize = 8),
        border = FALSE,
        gap = unit(c(1.5,0.5,0.5,0.5,0.5), "mm"),
        gp = gpar(col = "white",lty = "solid", lwd = 2.5)
        )

Bottom_column =  HeatmapAnnotation(
        Fraction = anno_barplot(fraction %>% select(`C>A`:`T>G`),bar_width = 0.7,
            gp = gpar(fill = c("#A2666F","#DCAE8D","#384A69","#D06C6B","#665775","#8E5C77"),
                    col = c("#A2666F","#DCAE8D","#384A69","#D06C6B","#665775","#8E5C77"), alpha = 1),
            height = unit(1.5, "cm"),border = FALSE),
           annotation_name_side = "left",
        annotation_name_gp = gpar(fontsize = 12)
)


lgd_list = list(

    Legend(labels = c("C>A","C>G", "C>T","T>C","T>A","T>G"),
        legend_gp = gpar(fill = c("#A2666F","#DCAE8D","#384A69","#D06C6B","#665775","#8E5C77")), 
        title = "Substitution", nr = 6,
        grid_height = unit(1, "mm"),grid_width = unit(4, "mm"))
)

pdf("test.pdf", width = 7, height = 7)
ht<- oncoPrint(df1,
    alter_fun = alter_fun, 
    col = col,
    width = unit(6, "cm"),
    height = unit(7, "cm"),
    show_pct = TRUE,
    show_column_names = TRUE,
    show_row_names = TRUE,
    row_names_side = "left",
    pct_side = "right",
    column_split = df0$phase,
    # row_labels = rownames(df1),
    bottom_annotation = Bottom_column,
    top_annotation = Top_column
    )
draw(ht,annotation_legend_list = lgd_list
)
dev.off()

