
Div <- data.frame()
for (i in c("Pre","Post")){
  for (x in c("P01","P04","P07")){
  for (j in Sample$majorCluster %>% unique){
      
      test <- scRNA_sub_T@meta.data %>% 
        filter( status == i & donor == x & Tcelltype == j) %>%
        select(cdr3,frequency) %>% 
        group_by(cdr3) %>% summarise(frequency = n()) %>%
        arrange(desc(frequency)) %>% drop_na()

      shannon_entropy <- vegan::diversity(test$frequency , "shannon")
      Pielou_index <- shannon_entropy / log(vegan::specnumber(test$frequency))
      Clonality = 1 - Pielou_index

      Gini.Simpson <- vegan::diversity(test$frequency, "simpson")
      inverse.Simpson <- vegan::diversity(test$frequency, "invsimpson")
      Gini_coefficient <- Gini(test$frequency)

      diversity <- data.frame(status = i,
                           donor = x,
                           majorCluster = j,
                           shannon_entropy = shannon_entropy,
                           Gini_Simpson = Gini.Simpson,
                           inv_Simpson = inverse.Simpson,
                           Clonality = Clonality,
                           Gini_coefficient = Gini_coefficient)
      Div <- rbind(Div, diversity)
    }
  }
}
