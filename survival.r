
fit <- survfit(Surv(DFS, status1) ~ driver, data = df)
fit <- survfit(Surv(OS, status2) ~ driver, data = df)


ggsurv <- ggsurvplot(
    fit, 
    data =  df,
    size = 1,
    palette = c("#5BBCD6","#EB3323"),
    pval = TRUE, 
    risk.table = TRUE, 
    xlab = "Time in Months", 
    ylab = "Overall Survival",
    # ylab="Disease-Free Survival",
    # legend.labs = c("Low","High"), # Change legend labels
    legend.labs = c("Other\ngene","Driver\ntrunk"), # Change legend labels
    # risk.table.title = "Disease-Specific Survival Event",
    risk.table.title = "Overall Survival",
    ggtheme = theme_bw()
    )
