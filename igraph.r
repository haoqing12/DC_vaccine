library(tidyverse)
library(data.table)
library(igraph)
library(ggraph)
library(stringr)
library(reshape2)

df <- fread("/work/haoq/NEO_scRNA/GLIPH/output/LDP_pre_post_new.tsv") %>% 
    filter(pattern != "single" & Fisher_score < 0.05) %>%
    select(index,pattern, type,number_unique_cdr3, ulTcRb, TcRb,expansion_score,Sample,Freq)

Parttern <- df %>% ungroup() %>% select(pattern) %>% 
    distinct()


n = nrow(Parttern)
m = matrix(0,n,n,dimnames=list(Parttern$pattern,Parttern$pattern))

connect <- m %>% as.data.frame() %>% tibble::rownames_to_column() %>% rename("from" = 1) %>%
  tidyr::gather(key="to", value="value", -1)

for (i in df[duplicated(df$TcRb)] %>% pull(TcRb)){
    tmp <- df %>% filter(TcRb== i) %>% pull(pattern)
    connect$value[connect$from %in% tmp & connect$to %in% tmp] <- 1
}

# connect %>% filter(from %in% tmp & connect$to %in% tmp)
L <- connect %>% filter(value == 1) %>% pull(from) %>% unique
O <- setdiff(connect %>% pull(from) %>% unique, L)

g <- graph_from_data_frame(
  connect %>% filter(value >0), 
  directed = TRUE, 
  vertices = df %>% group_by(pattern,type,number_unique_cdr3) %>% 
                summarise(Sample = paste0(unique(Sample) %>% sort(), collapse = ","), Freq = sum(Freq)) %>%
                mutate(tag = ifelse(pattern %in% L, pattern, NA)) %>% distinct
)

V(g)$Sample <- factor(V(g)$Sample, 
    levels = c("ES0-07:Pre","ES0-07:Post","ES0-07:Post,ES0-07:Pre")
    )

mycolor <- wes_palette("Chevalier1", 3, type = "discrete")
mycolor <- sample(mycolor, length(mycolor))


### expanded TCR network
g <- graph_from_data_frame(
  connect %>% filter(value >0), 
  directed = TRUE, 
  vertices = df %>% group_by(pattern,type,number_unique_cdr3) %>% 
                summarise(Sample = paste0(unique(Sample) %>% sort(), collapse = ","), 
                 Freq = sum(Freq),
                  expansion_score = min(expansion_score)) %>%
                mutate(Clonal = ifelse(expansion_score < 0.05, "YES", "No")) %>% distinct
)

## IFN G positive network
IFNG_exp_cells <- fread("IFNG_exp_cells.tsv")

V1 <- df %>% 
    left_join(scRNA_sub_T@meta.data %>% 
        select(barcode, cdr3, frequency,clonality_class, Tcelltype) %>%
        separate(cdr3,into=c('TRA','TRB'),sep=';') %>%
        separate(TRB,into=c(NA,'TRB_CDR3.aa'),sep=':') %>%
        separate(TRA,into=c(NA,'TRA_CDR3.aa'),sep=':'), by = c("TcRb"="TRB_CDR3.aa")) %>%
    left_join(IFNG_exp_cells, by = c("barcode"="."))

V1$IFNG[is.na(V1$IFNG)] <- 0

V2 <- V1 %>% group_by(pattern,type,number_unique_cdr3) %>% 
    summarise(IFNG = max(IFNG))

V3 <- V2 %>% ungroup %>% 
  left_join(df %>% group_by(pattern,type,number_unique_cdr3) %>% 
                summarise(Freq = sum(Freq)), by = c("pattern","type","number_unique_cdr3"))

V3 %>% filter(pattern == "VAYR")

g <- graph_from_data_frame(
  connect %>% filter(value >0), 
  directed = TRUE, 
  vertices = V3
)

coord <- layout_nicely(g) 

p<-ggraph(g,layout = "manual", x = coord[,1], y = coord[,2] ) + 
# p<-ggraph(g,layout = "fr") + 
  geom_edge_link(edge_colour="black", 
      edge_alpha=0.9, edge_width=0.5) +
  geom_node_point(aes(
    # size=V(g)$Freq,
     fill=as.factor(Sample)
    # fill = as.factor(Clonal)
    # fill = as.factor(IFNG)
    ), 
    shape=21,color='black',alpha=1)+
  geom_node_label(
    aes(label = V(g)$tag),
        size = 2, 
        fill = NA,color = "black",
        label.padding = unit(0, "mm"),
        label.size = unit(0, "mm"),
        repel = FALSE
        )+
  scale_size_continuous(range=c(6)) +
  
  scale_fill_manual(values=c("#D3DDDC","#FDD262","#446455"))+
  theme_minimal() +
  guides(fill = 'none')+
  theme(
    panel.grid = element_blank(),
    axis.line = element_blank(),
    axis.ticks =element_blank(),
    axis.text =element_blank(),
    axis.title = element_blank(),
  )
