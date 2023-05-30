library("scpdata")
library("ggplot2")
library("patchwork")
library("scp")
library("limma")
source("utils.R")

## For KNN imputation
k <- 15
## Missing protein filtering
pNA <- 0.99

####---- Load data ----####

dat <- leduc2022()
dat <- filterNA(dat, i = "proteins_norm2", pNA = pNA)

####---- Batch correction ----####

dat <- correctBatch(dat, i = "proteins_norm2", 
                    name = "proteins_bc", biolCol = "SampleType",
                    batchCol1 = "lcbatch", batchCol2 = "Channel")

####---- Imputation ----####

i <- "proteins_bc"
## KNN by row
dat <- impute(dat, i = i, name = "knn_row", MARGIN = 1,
              method = "knn", k = k, rowmax = 1, colmax = 1, 
              maxp = Inf)
## KNN by col
dat <- impute(dat, i = i, name = "knn_col", MARGIN = 2, 
              method = "knn", k = k, rowmax = 1, colmax = 1, 
              maxp = Inf)

knnMethods <- c("sample-wise", "none", "variable-wise")
knnlist <- list(getWithColData(dat, "knn_col"),
                getWithColData(dat, i),
                getWithColData(dat, "knn_row"))
names(knnlist) <- knnMethods

####---- Correlations ----####

## Compute cell correlations
cellCordf <- lapply(names(knnlist), function(ll) {
    x <- knnlist[[ll]]
    out <- correlationTable(assay(x), MARGIN = 2, 
                            colData(x)[, "SampleType", drop = FALSE])
    out$type <- ll
    out
})
cellCordf <- do.call(rbind, cellCordf)
cellCordf$type <- factor(cellCordf$type, levels = knnMethods)
cellCordf$SampleType <- recode(cellCordf$SampleType,
                               `within Melanoma cell` = "Melanoma",
                               `within Monocyte` = "Monocyte")

## Compute protein correlations
## First generate clusters of proteins
## Bin by percentage missing
pNAbinned <- binProteinsByPercentMissing(knnlist$none, ngroups = 5)
npNAbinned <- table(pNAbinned)

protCordf <- lapply(names(knnlist), function(ll) {
    x <- knnlist[[ll]]
    out <- correlationTable(assay(x), MARGIN = 1,
                            data.frame(pNAbinned = pNAbinned,
                                       row.names = names(pNAbinned)))
    out$type <- ll
    out
})
protCordf <- do.call(rbind, protCordf)

####---- Plot correlations ----####

(cellCorPlByGroup <- cellCordf %>% 
     filter(!grepl("between", SampleType)) %>%
     ggplot() +
     aes(y = cor,
         fill = SampleType,
         x = type) +
     geom_violin() +
     stat_summary(fun = median, colour = "grey20") +
     facet_grid(SampleType ~., scales = "free_x") +
     labs(y = "Cell correlation", x = "") +
     theme_minimal())
(cellCorAllPl <- cellCordf %>% 
        ggplot() +
        aes(y = cor,
            x = type) +
        geom_violin(fill = "grey") +
        stat_summary(fun = median, colour = "grey20") +
        labs(y = "Cell correlation", x = "") +
        theme_minimal())

# Protein correlation are unnecessarily large, thus we subsample
subsamp <- sample(seq_len(nrow(protCordf)), nrow(protCordf) / 100)
(protCorAllPl <- protCordf[subsamp, ] %>% 
        mutate(type = factor(type, levels = knnMethods)) %>% 
        ggplot() +
        aes(y = cor,
            x = type) +
        geom_violin(fill = "grey80") +
        stat_summary(fun = median, colour = "grey20") +
        labs(y = "Protein correlation", x = "") +
        theme_minimal())

(protCorByClusterPl <- protCordf[subsamp, ] %>% 
        filter(grepl("within", pNAbinned)) %>%
        mutate(pNAbinned = sub("within ", "", pNAbinned),
               pNAbinned = paste0("%missing: ", pNAbinned, "\n n = ",
                                  npNAbinned[pNAbinned])) %>%
        mutate(type = factor(type, levels = knnMethods)) %>% 
        ggplot() +
        aes(y = cor,
            fill = pNAbinned,
            x = type) +
        geom_violin() +
        stat_summary(fun = median, colour = "grey20") +
        scale_fill_manual(values = colorRampPalette(c("#6b8a58", "beige"))(5)) +
        facet_grid(pNAbinned ~ ., scales = "free_x") +
        labs(y = "Protein correlation", x = "") +
        theme_minimal())

#### ---- Make Figure ----####

(fig <- cellCorAllPl +
    protCorAllPl +
    cellCorPlByGroup +
    guides(fill = "none") +
    scale_fill_manual(values = c("#FF5733", "#048ABF")) +
    protCorByClusterPl +
    guides(fill = "none") +
    plot_layout(design = "14\n24\n34\n34") +
    plot_annotation(tag_levels = list(c("a", "b", "c", "d"))))

if (!dir.exists("figs")) {
    dir.create("figs")
}
ggsave("figs/figure3.pdf", fig, width = 7, height = 7.5)
