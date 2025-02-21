dir <- file.path(paste(info$sample_outs, "count/sample_feature_bc_matrix", sep="/"))
names(dir) = info$origin

counts <- Read10X(data.dir = dir)
scRNA_raw <- CreateSeuratObject(counts, min.cells = 3, min.features = 200)
scRNA_raw@meta.data %>% head()
table(scRNA_raw@meta.data$orig.ident)

scRNA_raw[["percent.mt"]] <- PercentageFeatureSet(object = scRNA_raw, pattern = "^MT-")

VlnPlot(object = scRNA_raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

scRNA_filter <- subset(scRNA_raw, subset = nFeature_RNA > 200 & nFeature_RNA < 4500 & percent.mt < 20)

#=============VDJ=============
dir_clone <- file.path(paste(info$sample_outs, "vdj_t/", sep="/"))

gettcr <- function(input, sample){
    tcr <- read.csv(paste(input,"filtered_contig_annotations.csv", sep="")) %>%
        #过滤1
        filter(productive == "true" & high_confidence == "true" & full_length =="true") %>%
        mutate(barcode = paste0(sample,"_",barcode)) %>%
        #### If a cell had two or more qualified chains of the same type, only that chain with the highest UMI count was qualified and retained.
        group_by(barcode,chain) %>%
        slice(which.max(umis)) %>%
        ungroup() %>% 
        mutate(cdr3 = paste0(chain,":",cdr3)) %>%
        group_by(barcode) %>%
        summarise(
            v_gene = paste(v_gene, collapse = ";"),
            j_gene = paste(j_gene, collapse = ";"),
            cdr3 = paste(cdr3, collapse = ";"),
            chain = n(),
            raw_clonotype_id = paste(raw_clonotype_id, collapse = ";")
        ) %>% 
        filter(chain > 1) %>% 
        select(-chain)

    freq <- tcr %>% ungroup() %>%
        group_by(cdr3) %>% 
        summarise(frequency = n()) %>%
        arrange(desc(frequency)) %>%
        mutate(
            new_clonotype_id = paste0(sample,"_clonotype",rownames(.)),
            proportion = frequency/sum(frequency)
        )
    
    TCR <- tcr %>% left_join(freq, by = "cdr3") %>% as.data.frame()
    return(TCR)
}

rownames(TCR) <- TCR$barcode
scRNA_filter <- AddMetaData(scRNA_filter, metadata = TCR)

num_PC <- 30
scRNA_harmony <- NormalizeData(scRNA_filter) %>% FindVariableFeatures() %>% ScaleData() %>% RunPCA(verbose=FALSE)
scRNA_harmony@meta.data$orig.ident <- as.factor(scRNA_harmony@meta.data$orig.ident)
Idents(scRNA_harmony) <- scRNA_harmony@meta.data$orig.ident
library(harmony)
scRNA_harmony <- RunHarmony(scRNA_harmony, group.by.vars = "orig.ident")
scRNA_harmony <- RunUMAP(scRNA_harmony, reduction = "harmony", dims = 1:num_PC) %>%
                FindNeighbors(reduction = "harmony", dims = 1:num_PC) %>% 
                FindClusters(resolution = 1)

scRNA_harmony.markers <- FindAllMarkers(scRNA_harmony, only.pos = TRUE, min.pct = 0.1, 
        logfc.threshold = 0.25)